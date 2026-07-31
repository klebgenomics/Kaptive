r"""Surface antigen serotyping and locus matching engine.

The `kaptive.serotyping` sub-package implements surface polysaccharide locus
typing (such as Klebsiella K and O loci, Acinetobacter K and OC loci) by matching
assembly contigs against reference locus databases, scoring gene presence,
and evaluating locus integrity.

**Exports:**

- **Serotyper**: Main serotyping execution engine
  ([`Serotyper`][kaptive.serotyping.core.Serotyper]).
- **SerotypingProblem**: Problem definition pairing a genome assembly with a database
  ([`SerotypingProblem`][kaptive.serotyping.models.SerotypingProblem]).
- **SerotypingResult**: Comprehensive outcome of locus typing and gene scoring
  ([`SerotypingResult`][kaptive.serotyping.models.SerotypingResult]).
- **GeneState**: Enumeration of gene call integrity states
  ([`GeneState`][kaptive.serotyping.models.GeneState]).
- **GeneHits**: Container for mapped gene alignment records
  ([`GeneHits`][kaptive.serotyping.models.GeneHits]).
- **LocusPieces**: Container for assembly locus fragment matches
  ([`LocusPieces`][kaptive.serotyping.models.LocusPieces]).
- **ReportRow**: Base container for structured output report serialization
  ([`ReportRow`][kaptive.serotyping.io.ReportRow]).
- **KaptiveRow**: Standard Kaptive TSV report format row
  ([`KaptiveRow`][kaptive.serotyping.io.KaptiveRow]).
- **Pha4geRow**: PHA4GE-compliant tabular report format row
  ([`Pha4geRow`][kaptive.serotyping.io.Pha4geRow]).
"""

from .core import Serotyper
from .io import KaptiveRow, Pha4geRow, ReportRow
from .models import GeneHits, GeneState, LocusPieces, SerotypingProblem, SerotypingResult

__all__ = [
    "GeneState",
    "SerotypingProblem",
    "GeneHits",
    "LocusPieces",
    "SerotypingResult",
    "Serotyper",
    "ReportRow",
    "KaptiveRow",
    "Pha4geRow",
]

