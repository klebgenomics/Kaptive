import os
import json
import uuid
import mimetypes
import gzip
from typing import Any
from urllib.parse import quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class KaptiveWebClientError(Exception):
    pass


# Classes --------------------------------------------------------------------------------------------------------------
class KaptiveWebClient:
    def __init__(self, api_key: str, base_url: str = "http://127.0.0.1:8000"):
        self.api_key = api_key
        # Ensure base_url does not end with a slash, and append /api if not present
        self.base_url = base_url.rstrip("/")
        if not self.base_url.endswith("/api"):
            self.base_url += "/api"

    def _request(self, endpoint: str, method: str = "GET", data: bytes = None, headers: dict[str, str] = None) -> Any:
        url = f"{self.base_url}{endpoint}"
        req_headers = {
            "X-API-Key": self.api_key,
            "Accept": "application/json",
        }
        if headers:
            req_headers.update(headers)

        req = Request(url, data=data, method=method, headers=req_headers)
        
        try:
            with urlopen(req) as response:
                body = response.read()
                if response.info().get("Content-Encoding") == "gzip":
                    body = gzip.decompress(body)  # If the response is gzipped, decompress it
                if response.info().get_content_type() == "application/json":
                    return json.loads(body.decode("utf-8"))
                return body

        except HTTPError as e:
            try:
                err_data = json.loads(e.read().decode("utf-8"))
                detail = err_data.get("detail", str(e))
            except Exception:
                detail = str(e)
            raise KaptiveWebClientError(f"HTTP {e.code}: {detail}")

        except URLError as e:
            raise KaptiveWebClientError(f"Network error: Failed to connect to {self.base_url}. Ensure you have an active internet connection. ({e.reason})")

    @staticmethod
    def _build_multipart_form(files: list[str]) -> tuple[bytes, str]:
        """Builds a multipart/form-data body using standard libraries."""
        boundary = uuid.uuid4().hex
        body = bytearray()
        
        for file_path in files:
            filename = os.path.basename(file_path)
            mime_type, _ = mimetypes.guess_type(file_path)
            if not mime_type:
                mime_type = "application/octet-stream"
                
            body.extend(f"--{boundary}\r\n".encode("utf-8"))
            body.extend(f'Content-Disposition: form-data; name="files"; filename="{filename}"\r\n'.encode("utf-8"))
            body.extend(f"Content-Type: {mime_type}\r\n\r\n".encode("utf-8"))
            with open(file_path, "rb") as f:
                body.extend(f.read())
            body.extend(b"\r\n")
            
        body.extend(f"--{boundary}--\r\n".encode("utf-8"))
        content_type = f"multipart/form-data; boundary={boundary}"
        return bytes(body), content_type

    def submit_genomes(self, species: str, files: list[str]) -> str:
        """Submits genomes for serotyping and returns the run_id."""
        body, content_type = self._build_multipart_form(files)
        # Endpoint expects URL-encoded species name
        endpoint = f"/serotype/{quote(species)}"
        response = self._request(endpoint, method="POST", data=body, headers={"Content-Type": content_type})
        return response.get("run_id")

    def get_run(self, run_id: str) -> dict[str, Any]:
        """Fetches the status and results of a run."""
        endpoint = f"/serotype/runs/{run_id}"
        return self._request(endpoint, method="GET")

    def download_jsonl(self, genome_ids: list[str]) -> bytes:
        """Downloads gzipped jsonl results for given genome IDs."""
        endpoint = "/serotype/results/download/jsonl"
        data = json.dumps({"genome_ids": genome_ids}).encode("utf-8")
        return self._request(endpoint, method="POST", data=data, headers={"Content-Type": "application/json"})
