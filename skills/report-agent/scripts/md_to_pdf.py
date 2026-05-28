#!/usr/bin/env python3
"""
SCAI Report Generator — Markdown to PDF via Pandoc + Typst

Converts SCAI markdown reports to publication-ready PDFs.
Uses pandoc for MD → Typst conversion and typst for Typst → PDF compilation.

Requires:
  - pandoc >= 3.1 (with typst writer)
  - typst >= 0.12

Usage:
  python md_to_pdf.py report.md -o report.pdf
  python md_to_pdf.py report.md -o report.pdf --template template.typ --margin 1in
"""

import argparse
import subprocess
import sys
from pathlib import Path


def check_dependencies() -> None:
    """Verify pandoc and typst are available."""
    missing = []

    try:
        subprocess.run(["pandoc", "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        missing.append("pandoc (install: conda install -c conda-forge pandoc)")

    try:
        result = subprocess.run(["typst", "--version"], capture_output=True, text=True)
        # Check typst can compile (version output means it works)
        if result.returncode != 0:
            missing.append("typst (install: see https://github.com/typst/typst)")
    except FileNotFoundError:
        missing.append("typst (install: see https://github.com/typst/typst)")

    if missing:
        print("❌ Missing dependencies:", file=sys.stderr)
        for m in missing:
            print(f"   • {m}", file=sys.stderr)
        sys.exit(1)


def md_to_typst(md_path: Path, typ_path: Path, template: Path | None = None) -> None:
    """Convert Markdown to Typst using pandoc."""
    cmd = ["pandoc", str(md_path), "-o", str(typ_path), "--to", "typst"]

    if template and template.exists():
        # Note: pandoc's --template for typst requires the template
        # to be in pandoc's template format, not a raw .typ file.
        # For now, we use it as a reference include via stdout manipulation.
        print(f"  📄 Using Typst template: {template}", file=sys.stderr)

    print(f"  🔄 Converting MD → Typst: {md_path.name} → {typ_path.name}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"❌ Pandoc error:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    if result.stderr:
        # Pandoc warnings (not errors)
        for line in result.stderr.strip().split("\n"):
            if line:
                print(f"  ⚠️  {line}", file=sys.stderr)


def fix_typst_syntax(content: str) -> str:
    """Fix common pandoc→Typst compatibility issues for Typst >= 0.12."""
    fixes = [
        # Typst 0.12+ uses #line() instead of #horizontalrule
        ("#horizontalrule", "#line(length: 100%)"),
        # Fix table column alignment syntax if needed
    ]
    for old, new in fixes:
        content = content.replace(old, new)
    return content


def typst_to_pdf(typ_path: Path, pdf_path: Path, margin: str = "1in") -> None:
    """Compile Typst to PDF."""
    print(f"  🔄 Compiling Typst → PDF: {typ_path.name} → {pdf_path.name}", file=sys.stderr)

    # Inject margin and font settings into the .typ file header
    # Typst uses native syntax — no quotes around lengths
    header = f'#set page(margin: {margin})\n#set text(font: "Liberation Sans", size: 11pt)\n'
    content = typ_path.read_text()

    # Fix pandoc compatibility issues
    content = fix_typst_syntax(content)

    # Only inject header if not already present (from previous runs)
    if "#set page" not in content[:300]:
        content = header + "\n" + content

    typ_path.write_text(content)

    # Run typst from the .typ file's directory so image paths resolve
    cwd = typ_path.parent
    result = subprocess.run(
        ["typst", "compile", str(typ_path.absolute()), str(pdf_path.absolute())],
        capture_output=True, text=True,
        cwd=str(cwd),
    )

    if result.returncode != 0:
        print(f"❌ Typst error:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    size = pdf_path.stat().st_size / 1024
    print(f"  ✅ PDF generated: {pdf_path} ({size:.0f} KB)", file=sys.stderr)


def md_to_pdf(
    md_path: str | Path,
    pdf_path: str | Path | None = None,
    template: str | Path | None = None,
    margin: str = "1in",
    keep_typ: bool = False,
) -> Path:
    """Convert a SCAI markdown report to PDF.

    Args:
        md_path: Path to input .md file.
        pdf_path: Path to output .pdf file. Default: same name as .md.
        template: Optional Typst template path (for reference).
        margin: Page margin for PDF (default: 1in).
        keep_typ: Keep intermediate .typ file (default: False).

    Returns:
        Path to the generated PDF.
    """
    md_path = Path(md_path)
    if not md_path.exists():
        print(f"❌ File not found: {md_path}", file=sys.stderr)
        sys.exit(1)

    if pdf_path is None:
        pdf_path = md_path.with_suffix(".pdf")
    else:
        pdf_path = Path(pdf_path)

    # Put .typ file next to the markdown so image paths resolve correctly
    typ_path = md_path.with_suffix(".typ")

    print(f"📄 Generating PDF from: {md_path}", file=sys.stderr)

    check_dependencies()
    md_to_typst(md_path, typ_path, Path(template) if template else None)
    typst_to_pdf(typ_path, pdf_path, margin)

    if not keep_typ:
        typ_path.unlink(missing_ok=True)
        print(f"  🧹 Removed intermediate: {typ_path.name}", file=sys.stderr)

    return pdf_path


def main():
    parser = argparse.ArgumentParser(
        description="SCAI Report Generator — Markdown to PDF via Pandoc + Typst",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Path to input Markdown file (.md)")
    parser.add_argument("-o", "--output", help="Path to output PDF file (default: same name as input)")
    parser.add_argument("--template", help="Path to Typst template (optional)")
    parser.add_argument("--margin", default="1in", help="Page margin (default: 1in)")
    parser.add_argument("--keep-typ", action="store_true", help="Keep intermediate .typ file")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages")

    args = parser.parse_args()

    if args.quiet:
        # Redirect stderr to /dev/null for the function calls
        pass

    pdf_path = md_to_pdf(
        md_path=args.input,
        pdf_path=args.output,
        template=args.template,
        margin=args.margin,
        keep_typ=args.keep_typ,
    )

    print(f"\n📄 Report ready: {pdf_path}")


if __name__ == "__main__":
    main()
