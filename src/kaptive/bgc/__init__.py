r"""Biosynthetic Gene Cluster (BGC) annotation, prediction, and model training.

The `kaptive.bgc` sub-package provides end-to-end tooling for bacterial BGC analysis:
fast ORF gene prediction and feature annotation, machine learning inference of cluster
architectures, and model training routines.

**Exports:**

- **Annotator**: ORF prediction and Randstrobe database search engine
  ([`Annotator`][kaptive.bgc.annotate.Annotator]).
- **AnnotationResult**: Container storing predicted genes, database hits, and BED export logic
  ([`AnnotationResult`][kaptive.bgc.annotate.AnnotationResult]).
- **Genes**: Structure-of-Arrays container for gene prediction batches
  ([`Genes`][kaptive.bgc.annotate.Genes]).
- **ArchitecturalModel**: Machine learning classifier architecture for BGC typing
  ([`ArchitecturalModel`][kaptive.bgc.models.ArchitecturalModel]).
- **ArchitecturalPredictor**: Inference engine predicting BGC regions from annotations
  ([`ArchitecturalPredictor`][kaptive.bgc.predictor.ArchitecturalPredictor]).
- **BGCPredictions**: Container for BGC region boundaries and probability scores
  ([`BGCPredictions`][kaptive.bgc.predictor.BGCPredictions]).
- **ArchitecturalTrainer**: Trainer for fitting BGC architectural models on training data
  ([`ArchitecturalTrainer`][kaptive.bgc.trainer.ArchitecturalTrainer]).
"""

from .annotate import AnnotationResult, Annotator, Genes
from .models import ArchitecturalModel
from .predictor import ArchitecturalPredictor, BGCPredictions
from .trainer import ArchitecturalTrainer

__all__ = [
    "Annotator",
    "AnnotationResult",
    "Genes",
    "ArchitecturalModel",
    "ArchitecturalPredictor",
    "BGCPredictions",
    "ArchitecturalTrainer",
]
