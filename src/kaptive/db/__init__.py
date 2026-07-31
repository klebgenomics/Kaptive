r"""Kaptive reference database management and data models.

The `kaptive.db` sub-package handles downloading, loading, validating, and
managing Kaptive reference databases containing locus sequences, gene metadata,
and serotype phenotype definitions.

**Exports:**

- **Database**: Core reference database containing indexed locus sequences
  ([`Database`][kaptive.db.core.Database]).
- **DatabaseManager**: Remote database downloading and local directory management
  ([`DatabaseManager`][kaptive.db.manager.DatabaseManager]).
- **DatabaseMetadata**: Metadata container for database versioning and locus definitions
  ([`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata]).
- **DatabaseError**: Base exception for database formatting and loading failures
  ([`DatabaseError`][kaptive.db.models.DatabaseError]).
- **Phenotype**: Single locus phenotype definition
  ([`Phenotype`][kaptive.db.models.Phenotype]).
- **Phenotypes**: Collection container for locus phenotype mapping records
  ([`Phenotypes`][kaptive.db.models.Phenotypes]).
"""

from .core import Database
from .manager import DatabaseManager
from .models import DatabaseError, DatabaseMetadata, Phenotype, Phenotypes

__all__ = [
    "DatabaseError",
    "DatabaseMetadata",
    "Phenotype",
    "Phenotypes",
    "Database",
    "DatabaseManager",
]

