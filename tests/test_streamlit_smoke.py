"""Headless smoke tests for every standalone Streamlit page."""

from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app/streamlit_app.py"
PAGES = (
    "Block overview",
    "Locus tracks",
    "Window tracks",
    "Explore a region",
    "Filter & export",
    "About & limitations",
)


def _navigation(app: AppTest):
    return next(item for item in app.sidebar.selectbox if item.label == "Go to")


@pytest.mark.parametrize("page", PAGES)
def test_each_page_renders_without_an_exception(page: str) -> None:
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    _navigation(app).set_value(page).run(timeout=60)
    assert not app.exception
    assert not app.error


def test_static_matplotlib_mode_renders() -> None:
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    renderer = next(
        item for item in app.sidebar.toggle if item.label == "Interactive plots (Plotly)"
    )
    renderer.set_value(False).run(timeout=60)
    assert not app.exception


def test_two_sessions_do_not_share_filter_state() -> None:
    first = AppTest.from_file(str(APP), default_timeout=60).run()
    second = AppTest.from_file(str(APP), default_timeout=60).run()
    _navigation(first).set_value("Filter & export").run(timeout=60)
    next(item for item in first.button if item.label == "Apply filters").click().run(timeout=60)
    assert not first.exception
    assert "filter_result" in first.session_state
    assert "filter_result" not in second.session_state
