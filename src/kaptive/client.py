r"""HTTP client for interacting with the Kaptive Web API.

This module provides [`KaptiveWebClient`][kaptive.client.KaptiveWebClient] for submitting
genome assemblies to remote Kaptive services, tracking job progress, and retrieving
serotyping result files.
"""

import gzip
import json
import mimetypes
import os
import uuid
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class KaptiveWebClientError(Exception):
    r"""Exception raised for HTTP network errors or API protocol failures."""

    pass


# Classes --------------------------------------------------------------------------------------------------------------
class KaptiveWebClient:
    r"""Client interface for sending requests to Kaptive Web API endpoints.

    Handles authentication header injection, request encoding, decompression of gzipped
    responses, and error handling.

    Attributes:
        api_key (str): Secret authentication API key string.
        base_url (str): Base URL of API service including `/api` path prefix.
    """

    def __init__(self, api_key: str, base_url: str = "http://127.0.0.1:8000") -> None:
        r"""Initialize web client with API key and server URL.

        Args:
            api_key (str): Secret API key token for authenticating requests.
            base_url (str): Host URL of the target API service. Defaults to `"http://127.0.0.1:8000"`.
        """
        self.api_key = api_key
        # Ensure base_url does not end with a slash, and append /api if not present
        self.base_url = base_url.rstrip("/")
        if not self.base_url.endswith("/api"):
            self.base_url += "/api"

    def _request(
        self,
        endpoint: str,
        method: str = "GET",
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
    ) -> Any:
        r"""Send HTTP request to specified API endpoint and parse response.

        Args:
            endpoint (str): API path suffix (e.g., `"/serotype/runs/123"`).
            method (str): HTTP request method verb (`"GET"`, `"POST"`, etc.). Defaults to `"GET"`.
            data (bytes | None): Raw request payload bytes to transmit. Defaults to `None`.
            headers (dict[str, str] | None): Optional dictionary of request header overrides. Defaults to `None`.

        Returns:
            Any: Decoded JSON dictionary/list or raw response bytes.

        Raises:
            KaptiveWebClientError: If an HTTP error response or network connectivity error occurs.
        """
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
            raise KaptiveWebClientError(
                f"Network error: Failed to connect to {self.base_url}. "
                f"Ensure you have an active internet connection. ({e.reason})"
            )

    @staticmethod
    def _build_multipart_form(files: list[str]) -> tuple[bytes, str]:
        r"""Build a multipart/form-data payload body from a list of file paths.

        Args:
            files (list[str]): List of local file system paths to include in form body.

        Returns:
            tuple[bytes, str]: Tuple containing raw encoded body bytes and content-type boundary string.
        """
        boundary = uuid.uuid4().hex
        body = bytearray()

        for file_path in files:
            filename = os.path.basename(file_path)
            mime_type, _ = mimetypes.guess_type(file_path)
            if not mime_type:
                mime_type = "application/octet-stream"

            body.extend(f"--{boundary}\r\n".encode())
            body.extend(f'Content-Disposition: form-data; name="files"; filename="{filename}"\r\n'.encode())
            body.extend(f"Content-Type: {mime_type}\r\n\r\n".encode())
            with open(file_path, "rb") as f:
                body.extend(f.read())
            body.extend(b"\r\n")

        body.extend(f"--{boundary}--\r\n".encode())
        content_type = f"multipart/form-data; boundary={boundary}"
        return bytes(body), content_type

    def submit_genomes(self, species: str, files: list[str]) -> str:
        r"""Submit genome assembly files for remote serotyping.

        Args:
            species (str): Target species name identifier.
            files (list[str]): List of file paths to genome assembly files.

        Returns:
            str: Remote execution run ID token string.

        Raises:
            KaptiveWebClientError: If request submission fails.
        """
        body, content_type = self._build_multipart_form(files)
        # Endpoint expects URL-encoded species name
        endpoint = f"/serotype/{quote(species)}"
        response = self._request(endpoint, method="POST", data=body, headers={"Content-Type": content_type})
        return response.get("run_id")

    def get_run(self, run_id: str) -> dict[str, Any]:
        r"""Fetch execution status and result dictionary for a serotyping job run.

        Args:
            run_id (str): Unique job execution run ID identifier.

        Returns:
            dict[str, Any]: Run status information and serotyping results dictionary.

        Raises:
            KaptiveWebClientError: If network request or server query fails.
        """
        endpoint = f"/serotype/runs/{run_id}"
        return self._request(endpoint, method="GET")

    def download_jsonl(self, genome_ids: list[str]) -> bytes:
        r"""Download JSONL serotyping results payload for specified genome IDs.

        Args:
            genome_ids (list[str]): List of genome identifiers to retrieve results for.

        Returns:
            bytes: Decompressed or raw bytes payload of JSONL result file.

        Raises:
            KaptiveWebClientError: If network request fails.
        """
        endpoint = "/serotype/results/download/jsonl"
        data = json.dumps({"genome_ids": genome_ids}).encode("utf-8")
        return self._request(endpoint, method="POST", data=data, headers={"Content-Type": "application/json"})

