from io import IOBase
from sys import stdout
from typing import Union, IO, Optional
from json import loads as json_loads
from zipfile import ZipFile
from tempfile import TemporaryDirectory
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from shutil import copyfileobj
from concurrent.futures import ThreadPoolExecutor, as_completed


# Classes --------------------------------------------------------------------------------------------------------------
# class LiteralFile(type(Path())):
#     """
#     A Path wrapper that evaluates to False if the file is missing or empty.
#     Inherits from the concrete Path type (PosixPath/WindowsPath) to ensure correct instantiation.
#     """
#     _MIN_SIZE = 1
#
#     def __bool__(self):
#         try: return self.is_file() and self.stat().st_size >= self._MIN_SIZE
#         except OSError: return False


class Downloader:
    """
    High-performance downloader supporting parallel chunk retrieval.
    """
    # _MAX_WORKERS = (Resources.available_cpus or 1) * 4
    _MAX_WORKERS = 32
    _CHUNK_SIZE = 10 * 1024 * 1024

    def __enter__(self):
        self._pool = ThreadPoolExecutor(max_workers=self._MAX_WORKERS)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if hasattr(self, '_pool'):
            self._pool.shutdown(wait=True)
            del self._pool

    def fetch(self, url: Union[str, Request], dest: Union[str, Path] = None, data=None,
              encode_data: bool = False) -> Union[Path, bytes]:
        if data is not None and encode_data: data = urlencode(data).encode('utf-8')
        if not isinstance(url, Request):
            url = Request(url, headers={'User-Agent': 'Mozilla/5.0'}, data=data)

        # Fallback for POST or small files
        if url.data is not None:
            return self._single_thread_download(url, dest)

        try:
            size, accept_ranges, real_url = self._get_info(url)
        except Exception:
            return self._single_thread_download(url, dest)

        if accept_ranges and size and size > self._CHUNK_SIZE:
            return self._parallel_download(real_url, dest, size, url.headers)

        return self._single_thread_download(url, dest)

    @staticmethod
    def _get_info(req: Request):
        head_req = Request(req.full_url, headers=req.headers, method='HEAD')
        with urlopen(head_req) as response:
            size = int(response.headers.get('Content-Length', 0))
            accept_ranges = response.headers.get('Accept-Ranges', 'none') == 'bytes'
            return size, accept_ranges, response.geturl()

    @staticmethod
    def _single_thread_download(url, dest):
        with urlopen(url) as response:
            if dest:
                dest = Path(dest)
                with open(dest, 'wb') as f: copyfileobj(response, f)
                return dest
            return response.read()

    def _parallel_download(self, url_str, dest, size, headers):
        chunks = []
        for i in range(0, size, self._CHUNK_SIZE):
            start = i
            end = min(i + self._CHUNK_SIZE - 1, size - 1)
            chunks.append((start, end))

        if dest:
            dest = Path(dest)
            # Create sparse file
            with open(dest, 'wb') as f:
                f.truncate(size)

            def _worker(start, end):
                req = Request(url_str, headers=headers)
                req.add_header('Range', f'bytes={start}-{end}')
                with urlopen(req) as response:
                    data = response.read()
                with open(dest, 'r+b') as f:
                    f.seek(start)
                    f.write(data)
        else:
            buffer = bytearray(size)

            def _worker(start, end):
                req = Request(url_str, headers=headers)
                req.add_header('Range', f'bytes={start}-{end}')
                with urlopen(req) as response:
                    data = response.read()
                buffer[start:end + 1] = data

        # Use existing pool if available (Context Manager), else create one
        pool = getattr(self, '_pool', None)
        if pool:
            futures = [pool.submit(_worker, s, e) for s, e in chunks]
            for f in as_completed(futures): f.result()
        else:
            with ThreadPoolExecutor(max_workers=self._MAX_WORKERS) as pool:
                futures = [pool.submit(_worker, s, e) for s, e in chunks]
                for f in as_completed(futures): f.result()

        return dest if dest else bytes(buffer)


class GitRepo:
    """
    A context-aware class to represent a git repository.
    Handles temporary directories automatically via the 'with' statement.

    Examples:
        >>> with GitRepo('owner', 'repo') as repo:
        ...     print(repo.local_path)
    """
    _BASE_URL = 'https://api.github.com'
    __slots__ = ('owner', 'repo', 'branch', '_api_url', '_meta_cache', '_temp_dir_obj', 'local_path')

    def __init__(self, owner: str, repo: str, branch: str = 'main'):
        self.owner = owner
        self.repo = repo
        self.branch = branch
        self._api_url = f'{self._BASE_URL}/repos/{owner}/{repo}'
        # Internal state
        self._meta_cache = None
        self._temp_dir_obj = None  # Holds the TemporaryDirectory object
        self.local_path: Optional[Path] = None

    def __repr__(self):
        return f"<GitRepo {self.owner}/{self.repo} ({self.branch})>"

    def fetch(self, filename: str):
        url = f'{self._api_url}/contents/{filename}?ref={self.branch}'
        req = Request(url, headers={'Accept': 'application/vnd.github.v3.raw'})
        try:
            with urlopen(req) as response:
                return response.read()
        except Exception as e:
            raise e

    def fetch_many(self, filenames: list[str]) -> dict[str, bytes]:
        """Fetch multiple files from the repository concurrently."""
        results = {}
        with ThreadPoolExecutor(max_workers=min(len(filenames), 10)) as pool:
            future_to_name = {pool.submit(self.fetch, name): name for name in filenames}
            for f in as_completed(future_to_name):
                results[future_to_name[f]] = f.result()
        return results

    @property
    def metadata(self) -> dict:
        """
        Lazy-loaded property for repository metadata.
        Only downloads from API once per instance.
        """
        if self._meta_cache is None:
            req = Request(self._api_url, headers={'Accept': 'application/vnd.github.v3+json'})
            with urlopen(req) as response:
                self._meta_cache = json_loads(response.read())
        return self._meta_cache

    def clone(self) -> Path:
        """
        Downloads and extracts the repo to a temporary directory.
        Returns the Path to the actual source code (inside the extracted folder).
        """
        # If already cloned, just return the path
        if self.local_path and self.local_path.exists(): return self.local_path
        zip_url = f'{self._api_url}/zipball/{self.branch}'

        # Create the temp directory explicitly
        self._temp_dir_obj = TemporaryDirectory()
        base_temp_path = Path(self._temp_dir_obj.name)
        zip_path = base_temp_path / "repo.zip"
        with Downloader() as dl:
            dl.fetch(zip_url, zip_path)

        with ZipFile(zip_path, 'r') as zfile:
            zfile.extractall(base_temp_path)
            root_folder = zfile.namelist()[0].split('/')[0]
        zip_path.unlink()
        self.local_path = base_temp_path / root_folder
        return self.local_path

    def cleanup(self):
        if self._temp_dir_obj:
            self._temp_dir_obj.cleanup()
            self._temp_dir_obj = None
            self.local_path = None

    def __enter__(self):
        """Enter the runtime context related to this object."""
        self.clone()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the runtime context and cleanup temporary directory."""
        self.cleanup()


# Functions ------------------------------------------------------------------------------------------------------------
def check_file(file: str | Path, min_size: int = 1) -> Path:
    if isinstance(file, str):
        file = Path(file)

    if file.is_file() and file.stat().st_size >= min_size:
        return file

    raise FileNotFoundError(file)


def write_to_file_or_directory(path: Union[str, Path, IO], mode: str = 'at') -> Union[Path, IO]:
    """
    Writes to a file or creates a directory based on the provided path.

    If the path is '-' or 'stdout', it returns stdout.
    If the path has a suffix, it's treated as a file and opened for appending.
    If the path has no suffix, it's treated as a directory and created if it doesn't exist.

    :param path: The path to the file or directory.
    :param mode: The mode to open the file in if it's a file.
    :return: A file handle (IO) if a file is specified, or a Path object if a directory is specified.
    """
    if isinstance(path, IOBase):
        return path
    if path in {'-', 'stdout'}:  # If the path is '-', return stdout
        return stdout
    if not isinstance(path, Path):  # Coerce to Path object
        path = Path(path)
    if path.suffix:  # If the path has an extension, it's probably a file
        # NB: We can't use is_file or is_dir because it may not exist yet, `open()` will create or append
        return open(path, mode)  # Open the file
    else:
        path.mkdir(exist_ok=True, parents=True)  # Create the directory if it doesn't exist
    return path