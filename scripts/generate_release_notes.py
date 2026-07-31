import json
import urllib.request
from urllib.error import URLError, HTTPError
import sys

def generate_release_notes(output_path):
    repo = "klebgenomics/Kaptive"
    url = f"https://api.github.com/repos/{repo}/releases"
    
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "Zensical-Release-Notes-Fetcher"
    }

    req = urllib.request.Request(url, headers=headers)
    
    try:
        with urllib.request.urlopen(req) as response:
            releases = json.loads(response.read().decode("utf-8"))
            
            with open(output_path, 'wt') as doc:
                doc.write("""---
title: Release Notes
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/rocket
categories:
  - Development
---

""")
                for release in releases[:5]:
                    name = release.get("name") or release.get("tag_name")
                    date = release.get("published_at", "")[:10]  # Extracts YYYY-MM-DD
                    body = release.get("body", "")
                    
                    doc.write(f"# {name}\n")
                    doc.write(f"*Published on {date}*\n\n")
                    doc.write(f"{body}\n\n")
                    doc.write("---\n\n")

    except HTTPError as e:
        print(f"HTTP Error: {e.code} - {e.reason}")
        sys.exit(1)
    except URLError as e:
        print(f"Failed to reach the server: {e.reason}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

def main():
    generate_release_notes('docs/releases.md')

if __name__ == '__main__':
    main()
