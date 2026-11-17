"""Prove that the extracted runtime is independent of the analysis repository."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_standalone_package_imports_without_analysis_package(tmp_path: Path) -> None:
    command = (
        "import sys; "
        f"sys.path.insert(0, {str(ROOT / 'src')!r}); "
        "from importlib.abc import MetaPathFinder; "
        "exec('class Blocker(MetaPathFinder):\\n"
        " def find_spec(self, fullname, path=None, target=None):\\n"
        "  if fullname.startswith(\"aa_haplotypes\"): raise ImportError(fullname)'); "
        "sys.meta_path.insert(0, Blocker()); "
        "import aa_metrics_app; "
        "assert aa_metrics_app.__file__.startswith(" + repr(str(ROOT / "src")) + "); "
        "assert 'aa_haplotypes' not in sys.modules"
    )
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, "-c", command],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_runtime_files_have_no_source_repository_path() -> None:
    runtime_roots = (ROOT / "app", ROOT / "src", ROOT / ".streamlit")
    text_suffixes = {".py", ".json", ".toml", ".css"}
    offenders: list[str] = []
    for runtime_root in runtime_roots:
        for path in runtime_root.rglob("*"):
            if path.is_file() and path.suffix in text_suffixes:
                if "15-msp1" in path.read_text(encoding="utf-8"):
                    offenders.append(path.relative_to(ROOT).as_posix())
    assert offenders == []


def test_only_active_application_configuration_and_assets_are_present() -> None:
    config_names = {path.name for path in (ROOT / "app/config").iterdir() if path.is_file()}
    assert config_names == {
        "app-content.json",
        "data-source.json",
        "filter-presets.json",
        "logo-help.json",
        "metric-help.json",
        "mutation-help.json",
        "structure-help.json",
    }
    assert not (ROOT / ".git").exists()
    assert not (ROOT / "notebooks").exists()
    assert not (ROOT / "_data").exists()
    assert not (ROOT / "src/aa_haplotypes").exists()
    assert not list(ROOT.rglob("*.sqlite"))


def test_configured_images_resolve_inside_the_app() -> None:
    path = ROOT / "app/config/app-content.json"
    content = json.loads(path.read_text(encoding="utf-8"))
    image_paths = [content["app"].get("logo_path"), content["app"].get("logo_icon_path")]
    image_paths.extend(page.get("icon_path") for page in content["pages"].values())
    for relative in (item for item in image_paths if item):
        resolved = (path.parent / relative).resolve()
        assert resolved.is_relative_to((ROOT / "app").resolve())
        assert resolved.is_file()


def test_extraction_manifest_covers_and_verifies_release_files() -> None:
    manifest_path = ROOT / "docs/extraction-manifest.tsv"
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    recorded = {row["target_path"]: row for row in rows}
    actual: set[str] = set()
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if relative.as_posix() == "docs/extraction-manifest.tsv":
            continue
        if any(part in {"__pycache__", ".pytest_cache", ".ruff_cache"} for part in relative.parts):
            continue
        if path.suffix == ".pyc":
            continue
        actual.add(relative.as_posix())
    assert set(recorded) == actual
    for relative, row in recorded.items():
        path = ROOT / relative
        assert int(row["size_bytes"]) == path.stat().st_size
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
