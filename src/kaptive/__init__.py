r"""Kaptive: Surface antigen locus serotyping and analysis toolkit.

This top-level package exposes version information and core module namespaces
for surface antigen serotyping, database management, and analysis.

Attributes:
    __version__ (str): Installed package version or "unknown" if not installed.
"""

try:
    import importlib.metadata

    __version__ = importlib.metadata.version("kaptive")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
