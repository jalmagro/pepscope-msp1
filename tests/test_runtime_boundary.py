"""Prove that the standalone runtime is independent of the analysis repository."""

from __future__ import annotations

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


def test_streamlit_media_server_dependencies_are_compatible() -> None:
    """Prevent Cloud builds from resolving an incompatible Starlette release."""

    requirements = {
        line.strip()
        for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    assert "streamlit==1.58.0" in requirements
    assert "starlette==1.4.1" in requirements


def test_configured_images_resolve_inside_the_app() -> None:
    path = ROOT / "app/config/app-content.json"
    content = json.loads(path.read_text(encoding="utf-8"))
    image_paths = [content["app"].get("logo_path"), content["app"].get("logo_icon_path")]
    image_paths.extend(page.get("icon_path") for page in content["pages"].values())
    for relative in (item for item in image_paths if item):
        resolved = (path.parent / relative).resolve()
        assert resolved.is_relative_to((ROOT / "app").resolve())
        assert resolved.is_file()
