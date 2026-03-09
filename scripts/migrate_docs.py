#!/usr/bin/env python3

# Imports --------------------------------------------------------------------------------------------------------------
from typing import Match
from re import MULTILINE, compile as regex
from subprocess import run
from pathlib import Path
from tempfile import NamedTemporaryFile
from pprint import pprint


# Constants ------------------------------------------------------------------------------------------------------------
_ANCHOR_REGEX = regex(r'\.\.\s+_([^:]+):')
_REF_TEXT_REGEX = regex(r':ref:`([^<]+?)\s*<([^>]+?)>`')
_REF_LABEL_REGEX = regex(r':ref:`([^`<]+?)`')
_ADMONITION_REGEX = regex(r'^> (?:\*\*([A-Za-z]+)\*\*|\[\!([A-Za-z]+)\])\n(?:>\n)?((?:> .*(?:\n|$))*)', flags=MULTILINE)


# Functions ------------------------------------------------------------------------------------------------------------
def slugify(text: str) -> str:
    """Creates a URL-safe slug from a label."""
    return text.lower().strip().replace(' ', '-').replace('_', '-')


def replace_ref_with_text(match: Match[str], label_map: dict[str, tuple[str, str]]) -> str:
    text = match.group(1).strip()
    label = match.group(2).strip().lower()
    if label in label_map:
        target_md, slug = label_map[label]
        return f'`{text} <{target_md}#{slug}>`_'
    return match.group(0)


def replace_target(match: Match[str]) -> str:
    label = match.group(1).strip()
    slug = slugify(label)
    return f'.. raw:: html\n\n   <a id="{slug}"></a>\n\n'


def replace_ref_label_only(match: Match[str], label_map: dict[str, tuple[str, str]]) -> str:
    label = match.group(1).strip()
    label_lower = label.lower()
    if label_lower in label_map:
        target_md, slug = label_map[label_lower]
        return f'`{label} <{target_md}#{slug}>`_'
    return match.group(0)


def format_admonition(match: Match[str]) -> str:
    # Get the type from whichever regex group matched (e.g., "Note" or "NOTE")
    adm_type = (match.group(1) or match.group(2)).lower()
    body = match.group(3)

    # Clean up the blockquote arrows and apply 4-space indent for MkDocs
    cleaned_body = []
    for line in body.splitlines():
        if line.startswith('> '):
            cleaned_body.append('    ' + line[2:])
        elif line.startswith('>'):
            cleaned_body.append('    ' + line[1:])
        else:
            cleaned_body.append('    ' + line)

    return f'!!! {adm_type}\n' + '\n'.join(cleaned_body) + '\n\n'


def build_label_map(rst_files: list[Path]) -> dict[str, tuple[str, str]]:
    label_map = {}
    for rst_file in rst_files:
        md_name = rst_file.with_suffix('.md').name
        for match in _ANCHOR_REGEX.finditer(rst_file.read_text()):
            label = match.group(1).strip()
            label_map[label.lower()] = (md_name, slugify(label))
    print(f"Found {len(label_map)} cross-reference targets.")
    return label_map


def convert(rst_file: Path, label_map: dict[str, tuple[str, str]]) -> Path:
    md_file = rst_file.with_suffix('.md')

    # Note: NamedTemporaryFile opens in binary mode by default, so we need mode='w'
    with NamedTemporaryFile(mode='w', encoding='utf-8', delete=False, suffix='.rst') as tmp_rst:
        tmp_rst_path = Path(tmp_rst.name)

        # PASS 2: Pre-process RST links
        content = rst_file.read_text(encoding='utf-8')

        # Fix 1: Inject raw HTML anchors
        content = _ANCHOR_REGEX.sub(replace_target, content)

        # Fix 2: Convert :ref:`Text <Label>` using lambda to pass label_map
        content = _REF_TEXT_REGEX.sub(lambda m: replace_ref_with_text(m, label_map), content)

        # Fix 3: Convert :ref:`Label` using lambda to pass label_map
        content = _REF_LABEL_REGEX.sub(lambda m: replace_ref_label_only(m, label_map), content)

        # Write temp file and run Pandoc
        tmp_rst.write(content)
        tmp_rst.flush()  # Ensure it's written to disk before Pandoc reads it

        print(f"Converting {rst_file} to {md_file}...")
        run(['pandoc', '-s', str(tmp_rst_path), '-t', 'gfm', '-o', str(md_file)])

    # Clean up the temp file (NamedTemporaryFile delete=False is safer across OS platforms)
    tmp_rst_path.unlink()

    # PASS 3: Post-process Markdown for Admonitions
    md_content = md_file.read_text(encoding='utf-8')
    md_content = _ADMONITION_REGEX.sub(format_admonition, md_content)
    md_file.write_text(md_content, encoding='utf-8')

    return md_file

def main():
    base_dir = Path('.')
    rst_files = list(base_dir.rglob('*.rst'))
    print(f"Found {len(rst_files)} rst files.")
    # PASS 1: Build a map of all cross-reference labels across all files
    label_map = build_label_map(rst_files)
    pprint(label_map)
    for rst_file in rst_files:
        convert(rst_file, label_map)

    print("\nConversion complete! Cross-references and admonitions have been preserved.")


# Entry point ----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()