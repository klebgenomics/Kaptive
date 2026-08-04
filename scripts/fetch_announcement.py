# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

import json
import os
import urllib.request

# Configuration
REPO_OWNER = "klebgenomics"
REPO_NAME = "Kaptive"
CATEGORY_NAME = "Announcements"

# GITHUB_TOKEN is automatically provided in GitHub Actions
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")


def main() -> None:
    query = f"""
    query {{
      repository(owner: "{REPO_OWNER}", name: "{REPO_NAME}") {{
        discussions(first: 20, orderBy: {{field: CREATED_AT, direction: DESC}}) {{
          nodes {{
            title
            url
            category {{
              name
            }}
          }}
        }}
      }}
    }}
    """

    if GITHUB_TOKEN:
        req = urllib.request.Request("https://api.github.com/graphql", method="POST")
        req.add_header("Authorization", f"Bearer {GITHUB_TOKEN}")
        req.add_header("Content-Type", "application/json")

        try:
            with urllib.request.urlopen(req, data=json.dumps({"query": query}).encode("utf-8")) as response:
                data = json.loads(response.read())

                # Filter the discussions to find the latest one in the target category
                discussions = data["data"]["repository"]["discussions"]["nodes"]
                announcement = next((d for d in discussions if d["category"]["name"] == CATEGORY_NAME), None)

                if announcement:
                    # Format however you want it to look in the banner!
                    html = f'<strong>📢 Latest Update:</strong> <a href="{announcement["url"]}">{announcement["title"]}</a>'

                    # Write to the file that main.html includes
                    with open("overrides/announcement.html", "w") as f:
                        f.write(html)

        except Exception as e:
            print(f"Failed to fetch announcement: {e}")


if __name__ == "__main__":
    main()
