"""Repository documentation: every local figure in README and gallery must exist."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def test_documentation_images_exist():
    for page in (ROOT / "README.md", ROOT / "docs" / "FIGURE_GALLERY.md"):
        for path in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", page.read_text(encoding="utf-8")):
            if path.startswith(("https://", "http://")):
                continue
            assert (page.parent / path).is_file(), (page, path)

def test_gallery_has_all_generated_pngs():
    contents = (ROOT / "docs" / "FIGURE_GALLERY.md").read_text(encoding="utf-8")
    for folder in ("figures", "validation_figures", "extended_figures", "thermo_evidence/figures"):
        for png in (ROOT / folder).glob("*.png"):
            assert png.name in contents
