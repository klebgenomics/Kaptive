# /// script
# requires-python = ">=3.11"
# dependencies = [deepl]
# ///

import os
import shutil
import re
from pathlib import Path
import deepl

LANGUAGES = {"ES": "es", "FR": "fr", "ZH": "zh", "JA": "ja"}


def translate_markdown(content: str, translator: deepl.Translator, target_lang: str) -> str:
    # Separate frontmatter
    frontmatter = ""
    body = content
    match = re.match(r"^(---\n.*?\n---\n)(.*)", content, re.DOTALL)
    if match:
        frontmatter = match.group(1)
        body = match.group(2)

        # Translate title and description in frontmatter
        for key in ["title", "description"]:
            key_match = re.search(rf"^({key}:\s*)(.*)$", frontmatter, re.MULTILINE)
            if key_match:
                prefix = key_match.group(1)
                value = key_match.group(2).strip("'\"")
                translated_val = translator.translate_text(value, target_lang=target_lang).text
                frontmatter = frontmatter.replace(key_match.group(0), f"{prefix}{translated_val}")

    if not body.strip():
        return frontmatter + body

    # We use tag_handling="xml" to prevent Deepl from messing up markdown as much, though it's experimental for markdown
    # Actually, plain text translation usually works well for Markdown in DeepL
    translated_body = translator.translate_text(body, target_lang=target_lang).text
    return frontmatter + translated_body


def main():
    api_key = os.environ.get("DEEPL_API_KEY")
    if not api_key:
        print("Error: DEEPL_API_KEY environment variable is not set.")
        print("Run with: DEEPL_API_KEY='your-key' just translate")
        exit(1)

    translator = deepl.Translator(api_key)

    docs_dir = Path("docs")

    for lang_code, lang_dir in LANGUAGES.items():
        print(f"Translating to {lang_code}...")
        target_dir = Path(f"docs_{lang_dir}")

        # Copy everything
        if target_dir.exists():
            shutil.rmtree(target_dir)
        shutil.copytree(docs_dir, target_dir)

        # Translate all md files
        for md_file in target_dir.rglob("*.md"):
            # Skip API reference files since they rely on mkdocstrings and Python source
            if "reference" in md_file.parts:
                continue

            with open(md_file, "r", encoding="utf-8") as f:
                content = f.read()

            translated_content = translate_markdown(content, translator, lang_code)

            with open(md_file, "w", encoding="utf-8") as f:
                f.write(translated_content)

        print(f"Finished translating {lang_code}")


if __name__ == "__main__":
    main()
