"""Validated, application-owned page and section copy."""

from __future__ import annotations

import json
import string
from dataclasses import dataclass
from html import escape
from pathlib import Path

APP_CONTENT_SCHEMA_VERSION = "1.3.0"
REQUIRED_PAGES = (
    "Block overview",
    "Locus tracks",
    "Window tracks",
    "Explore a region",
    "Filter & export",
    "About & limitations",
)
REQUIRED_SECTIONS = (
    "block_overview",
    "representative_block_sequences",
    "locus_tracks",
    "window_tracks",
    "peptide_sequence_inspector",
    "region_explorer",
    "matching_peptide_inspector",
    "mutation_explorer",
    "residue_logo",
    "structure_context",
    "filter_explorer",
    "filter_diagnostics",
    "conservation_target",
    "about",
    "metric_guide",
)
TITLE_TEMPLATE_FIELDS = frozenset({"app_title", "page_title"})


@dataclass(frozen=True)
class PageContent:
    """Title, introduction, and accessible local image for one navigation page."""

    title: str
    lead: str
    icon_path: Path
    icon_alt: str


@dataclass(frozen=True)
class SectionContent:
    """Heading and lead paragraph for one primary application section."""

    heading: str
    intro: str


@dataclass(frozen=True)
class AppContent:
    """Complete validated application copy catalog."""

    kicker: str
    title: str
    title_template: str
    sidebar_intro: str
    current_configuration_label: str
    favicon: str
    favicon_path: Path | None
    logo_path: Path | None
    logo_icon_path: Path | None
    pages: dict[str, PageContent]
    sections: dict[str, SectionContent]

    @property
    def page_icon(self) -> str | Path:
        """Return the configured image favicon or its symbol fallback."""

        return self.favicon_path or self.favicon

    def page(self, navigation: str) -> PageContent:
        """Return configured content for one navigation label."""

        try:
            return self.pages[navigation]
        except KeyError as error:
            raise ValueError(f"No application content for page {navigation!r}") from error

    def section(self, key: str) -> SectionContent:
        """Return configured content for one section key."""

        try:
            return self.sections[key]
        except KeyError as error:
            raise ValueError(f"No application content for section {key!r}") from error

    def page_title(self, navigation: str) -> str:
        """Return the active page's complete contextual title."""

        page = self.page(navigation)
        return self.title_template.format(
            app_title=self.title,
            page_title=page.title,
        )


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Return a mapping while rejecting duplicate JSON object keys."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Application content JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def _text(mapping: dict[str, object], field: str, context: str) -> str:
    value = mapping.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Application content {context}.{field} must be non-empty text")
    return value.strip()


def _text_with_default(
    mapping: dict[str, object], field: str, context: str, *, default: str
) -> str:
    """Return optional non-empty text, falling back only when the key is absent."""

    if field not in mapping:
        return default
    return _text(mapping, field, context)


def _optional_asset_path(
    mapping: dict[str, object],
    field: str,
    context: str,
    *,
    relative_to: Path,
) -> Path | None:
    """Validate an optional local image path and return its absolute path."""

    value = mapping.get(field)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            f"Application content {context}.{field} must be null or a non-empty path"
        )
    candidate = Path(value.strip())
    if not candidate.is_absolute():
        candidate = relative_to / candidate
    candidate = candidate.resolve()
    if not candidate.is_file():
        raise ValueError(
            f"Application content {context}.{field} does not identify an existing file: "
            f"{value!r}"
        )
    return candidate


def _page_icon_path(
    mapping: dict[str, object],
    field: str,
    context: str,
    *,
    content_path: Path,
) -> Path:
    """Return one validated PNG path contained within the application directory."""

    value = _text(mapping, field, context)
    relative = Path(value)
    if relative.is_absolute():
        raise ValueError(f"Application content {context}.{field} must be a relative path")
    app_root = (
        content_path.parent.parent
        if content_path.parent.name == "config"
        else content_path.parent
    )
    app_root = app_root.resolve()
    candidate = (content_path.parent / relative).resolve()
    try:
        candidate.relative_to(app_root)
    except ValueError as error:
        raise ValueError(
            f"Application content {context}.{field} must stay within {app_root}"
        ) from error
    if candidate.suffix.lower() != ".png":
        raise ValueError(f"Application content {context}.{field} must identify a PNG image")
    if not candidate.is_file():
        raise ValueError(
            f"Application content {context}.{field} does not identify an existing file: "
            f"{value!r}"
        )
    return candidate


def _exact_keys(mapping: dict[str, object], expected: tuple[str, ...], context: str) -> None:
    missing = sorted(set(expected).difference(mapping))
    unknown = sorted(set(mapping).difference(expected))
    if missing or unknown:
        raise ValueError(
            f"Application content {context} keys do not match the application. "
            f"Missing: {missing}; unknown: {unknown}"
        )


def load_app_content(path: str | Path) -> AppContent:
    """Load and validate editable application copy from JSON."""

    content_path = Path(path).resolve()
    payload = json.loads(
        content_path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
    )
    if not isinstance(payload, dict):
        raise ValueError("Application content JSON root must be an object")
    if payload.get("schema_version") != APP_CONTENT_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported application content schema version: "
            f"{payload.get('schema_version')!r}; expected {APP_CONTENT_SCHEMA_VERSION!r}"
        )
    app = payload.get("app")
    pages = payload.get("pages")
    sections = payload.get("sections")
    if not isinstance(app, dict) or not isinstance(pages, dict) or not isinstance(sections, dict):
        raise ValueError("Application content JSON requires app, pages, and sections objects")

    _exact_keys(pages, REQUIRED_PAGES, "pages")
    _exact_keys(sections, REQUIRED_SECTIONS, "sections")
    title_template = _text(app, "title_template", "app")
    fields = {
        field_name
        for _, field_name, _, _ in string.Formatter().parse(title_template)
        if field_name is not None
    }
    if fields != TITLE_TEMPLATE_FIELDS:
        raise ValueError(
            "Application title_template must use exactly {app_title} and {page_title}"
        )

    page_content: dict[str, PageContent] = {}
    for key in REQUIRED_PAGES:
        value = pages[key]
        if not isinstance(value, dict):
            raise ValueError(f"Application content page {key!r} must be an object")
        title = _text(value, "title", f"pages.{key}")
        lead = _text(value, "lead", f"pages.{key}")
        icon_path = _page_icon_path(
            value,
            "icon_path",
            f"pages.{key}",
            content_path=content_path,
        )
        icon_alt = _text(value, "icon_alt", f"pages.{key}")
        page_content[key] = PageContent(
            title=title,
            lead=lead,
            icon_path=icon_path,
            icon_alt=icon_alt,
        )

    section_content: dict[str, SectionContent] = {}
    for key in REQUIRED_SECTIONS:
        value = sections[key]
        if not isinstance(value, dict):
            raise ValueError(f"Application content section {key!r} must be an object")
        section_content[key] = SectionContent(
            heading=_text(value, "heading", f"sections.{key}"),
            intro=_text(value, "intro", f"sections.{key}"),
        )
    return AppContent(
        kicker=_text(app, "kicker", "app"),
        title=_text(app, "title", "app"),
        title_template=title_template,
        sidebar_intro=_text(app, "sidebar_intro", "app"),
        current_configuration_label=_text(app, "current_configuration_label", "app"),
        favicon=_text_with_default(app, "favicon", "app", default="🧬"),
        favicon_path=_optional_asset_path(
            app,
            "favicon_path",
            "app",
            relative_to=content_path.parent,
        ),
        logo_path=_optional_asset_path(
            app,
            "logo_path",
            "app",
            relative_to=content_path.parent,
        ),
        logo_icon_path=_optional_asset_path(
            app,
            "logo_icon_path",
            "app",
            relative_to=content_path.parent,
        ),
        pages=page_content,
        sections=section_content,
    )


def page_header_html(content: AppContent, navigation: str) -> str:
    """Return escaped semantic markup for the active page header."""

    return "".join(
        (
            page_kicker_html(content),
            page_title_html(content, navigation),
            page_lead_html(content, navigation),
        )
    )


def page_kicker_html(content: AppContent) -> str:
    """Return escaped markup for the small application kicker."""

    return f'<div class="app-kicker">{escape(content.kicker)}</div>'


def page_title_html(content: AppContent, navigation: str) -> str:
    """Return escaped markup for the active page title."""

    return f'<h1 class="app-title">{escape(content.page_title(navigation))}</h1>'


def page_lead_html(content: AppContent, navigation: str) -> str:
    """Return escaped markup for the active page's introduction."""

    return f'<p class="app-lede">{escape(content.page(navigation).lead)}</p>'


def section_intro_html(section: SectionContent) -> str:
    """Return escaped markup for a section's lead paragraph."""

    return f"<div class='section-intro'>{escape(section.intro)}</div>"
