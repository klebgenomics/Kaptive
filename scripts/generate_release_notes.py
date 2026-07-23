import json
import urllib.request
from urllib.error import URLError, HTTPError

# Replace with your GitHub repository details
repo = "owner/repo-name"
url = f"https://api.github.com/repos/{repo}/releases"

# GitHub API requires a User-Agent header
headers = {
    "Accept": "application/vnd.github+json",
    "User-Agent": "Zensical-Release-Notes-Fetcher"
}

# If your repo is private, uncomment the line below and add your token
# headers["Authorization"] = "Bearer YOUR_PERSONAL_ACCESS_TOKEN"

req = urllib.request.Request(url, headers=headers)

try:
    with urllib.request.urlopen(req) as response:
        # Read and parse the JSON payload
        releases = json.loads(response.read().decode("utf-8"))
        
        # Formats and prints the 5 most recent releases
        for release in releases[:5]:
            name = release.get("name") or release.get("tag_name")
            date = release.get("published_at", "")[:10]  # Extracts YYYY-MM-DD
            body = release.get("body", "")
            
            print(f"## {name}")
            print(f"*Published on {date}*\n")
            print(f"{body}\n")
            print("---")

except HTTPError as e:
    print(f"HTTP Error: {e.code} - {e.reason}")
except URLError as e:
    print(f"Failed to reach the server: {e.reason}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
