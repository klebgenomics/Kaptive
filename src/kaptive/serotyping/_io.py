# from warnings import warn
# from dataclasses import dataclass
# from io import TextIOBase
# from os import fstat
# from pathlib import Path
# from typing import Union, IO, Self
# from abc import ABC, abstractmethod
# from json import loads as json_loads, dumps as json_dumps
#
# # from kaptive.serotyping import GeneResult
# from kaptive.db import Database
#
#
# # Exceptions and warnings ----------------------------------------------------------------------------------------------
#
#
# # Classes -------------------------------------------------------------------------------------------------------------
# @dataclass
# class GenotypingSpec:
#     sample: str
#     genotyping_method: str
#     genotyping_schema_taxon: str
#     genotyping_database_name: str
#     genotyping_database_version: str
#     genotyping_schema_name: str
#     genotyping_software_name: str
#     genotyping_software_version: str
#     genotype: str
#     genotype_confidence_value: str
#     genotype_predicted_phenotype: str
#
#
# class GenotypingFormatter: pass
#
#
# class TypingResultFormatter:
#     _ASSEMBLY_HEADER = ('Assembly\tBest match locus\tBest match type\tMatch confidence\tProblems\tIdentity\tCoverage\t'
#                         'Length discrepancy\tExpected genes in locus\tExpected genes in locus, details\t'
#                         'Missing expected genes\tOther genes in locus\tOther genes in locus, details\t'
#                         'Expected genes outside locus\tExpected genes outside locus, details\t'
#                         'Other genes outside locus\tOther genes outside locus, details\t'
#                         'Truncated genes, details\tExtra genes, details\n')
#     _SCORES_HEADER = 'Assembly\tLocus\tAS\tmlen\tblen\tq_len\tgenes_found\tgenes_expected\n'
#
#     @classmethod
#     def from_dict(cls, d: dict, db: Database) -> TypingResult:
#         if not (best_match := db.loci.get(d['best_match'])):
#             raise TypingResultError(f"Best match {d['best_match']} not found in database")
#         self = TypingResult(sample_name=d['sample_name'], db=db, best_match=best_match,
#                             missing_genes=d['missing_genes'])
#         # Set the cached properties
#         self._percent_identity = float(d['percent_identity'])
#         self._percent_coverage = float(d['percent_coverage'])
#         self._phenotype = d['phenotype']
#         self._problems = d['problems']
#         self._confidence = d['confidence']
#         # Add the pieces and create the gene results
#         self.pieces = [LocusPiece.from_dict(i, result=self) for i in d['pieces']]
#         pieces = {i.__repr__(): i for i in self.pieces}
#         gene_results = {}  # This was previously a dict comp, but we need to check the gene is in the database, see #31
#         for r in chain(d['expected_genes_inside_locus'], d['unexpected_genes_inside_locus'],
#                        d['expected_genes_outside_locus'], d['unexpected_genes_outside_locus'],
#                        d['extra_genes']):
#             if not (gene := db.genes.get(r['gene'])) and not (gene := db.extra_genes.get(r['gene'])):
#                 raise TypingResultError(f"Gene {r['gene']} not found in database")
#             x = GeneResult.from_dict(r, result=self, piece=pieces.get(r['piece']), gene=gene)
#             gene_results[x.__repr__()] = x
#
#         for gene_result in gene_results.values():
#             self.add_gene_result(gene_result)
#         return self
#
#     @classmethod
#     def from_jsonl(cls, line: str, db: Database) -> Self:
#         try:
#             return TypingResult.from_dict(json_loads(line), db)
#         except Exception as e:
#             raise RuntimeError(f"Error parsing JSON line") from e
#
#     def format(self, format_spec) -> Union[str, dict]:
#         if format_spec == 'tsv':
#             return '\t'.join(
#                 [
#                     self.sample_name,
#                     self.best_match.name,
#                     self.phenotype,
#                     self.confidence,
#                     self.problems,
#                     f"{self.percent_identity:.2f}%",
#                     f"{self.percent_coverage:.2f}%",
#                     f"{self.__len__() - len(self.best_match)} bp" if len(self.pieces) == 1 else 'n/a',
#                     f"{(n_inside := len({i.gene.name for i in self.expected_genes_inside_locus}))} / "
#                     f"{(n_expected := len(self.best_match.seqs))} ({100 * n_inside / n_expected:.2f}%)",
#                     ';'.join(map(str, self.expected_genes_inside_locus)),
#                     ';'.join(self.missing_genes),
#                     f"{len(self.unexpected_genes_inside_locus)}",
#                     ';'.join(map(str, self.unexpected_genes_inside_locus)),
#                     f"{(n_outside := len({i.gene.name for i in self.expected_genes_outside_locus}))} / {n_expected} ({100 * n_outside / n_expected:.2f}%)",
#                     ';'.join(map(str, self.expected_genes_outside_locus)),
#                     f"{len(self.unexpected_genes_outside_locus)}",
#                     ';'.join(map(str, self.unexpected_genes_outside_locus)),
#                     ';'.join(map(str, filter(lambda z: z.phenotype == "truncated", self))),
#                     ';'.join(map(str, self.extra_genes))
#                 ]
#             ) + "\n"
#         if format_spec == 'fna':  # Return the nucleotide sequence of the locus
#             return "".join([i.format(format_spec) for i in self.pieces])
#         if format_spec in {'faa', 'ffn'}:  # Return the protein or nucleotide sequence of gene results
#             return "".join([i.format(format_spec) for i in self])
#         if format_spec == 'json':
#             return json_dumps(
#                 {
#                     'sample_name': self.sample_name, 'best_match': self.best_match.name, 'confidence': self.confidence,
#                     'phenotype': self.phenotype, 'problems': self.problems,
#                     'percent_identity': str(self.percent_identity),
#                     'percent_coverage': str(self.percent_coverage), 'missing_genes': self.missing_genes
#                 } | {
#                     attr: [i.format(format_spec) for i in getattr(self, attr)] for attr in {
#                         'pieces', 'expected_genes_inside_locus', 'unexpected_genes_inside_locus',
#                         'expected_genes_outside_locus', 'unexpected_genes_outside_locus', 'extra_genes'
#                     }
#                 }) + "\n"
#         raise ValueError(f"Unknown format specifier {format_spec}")
#
#     def write(self,
#               tsv: IO[str] = None,
#               json: IO[str] = None,
#               fna: Union[str, Path, IO[str]] = None,
#               ffn: Union[str, Path, IO[str]] = None,
#               faa: Union[str, Path, IO[str]] = None,
#               ):
#         """Write the serotyping result to files or file handles."""
#         [f.write(self.format(fmt)) for f, fmt in [(tsv, 'tsv'), (json, 'json')] if isinstance(f, TextIOBase)]
#         for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
#             if f:
#                 if isinstance(f, TextIOBase):
#                     f.write(self.format(fmt))
#                 elif isinstance(f, (Path, str)):
#                     Path(f / f'{self.sample_name}_kaptive_results.{fmt}').write_text(self.format(fmt))
#
#     def write_headers(tsv: TextIO = None, no_header: bool = False, scores: bool = False) -> int:
#         """Write appropriate header to a file handle."""
#         if tsv and not no_header and (tsv.name == '<stdout>' or fstat(tsv.fileno()).st_size == 0):
#             return tsv.write(_SCORES_HEADER if scores else _ASSEMBLY_HEADER)
#
#
# class LocusPieceFormatter:
#     @classmethod
#     def from_dict(cls, d: dict, **kwargs) -> LocusPiece:
#         return cls(id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'],
#                    sequence=Seq(d['sequence']), **kwargs)
#
#     def format(self, format_spec) -> Union[str, dict]:
#         if format_spec == 'fna':
#             return f">{self.result.genome}|{self.id}:{self.start}-{self.end}{self.strand}\n{self.sequence}\n"
#         if format_spec == 'json':
#             return {'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
#                     'sequence': str(self.sequence)}
#         raise ValueError(f"Unknown format specifier {format_spec}")
#
#
# class GeneResultFormatter:
#     def __str__(self) -> str:
#         s = f'{self.gene.name},{self.percent_identity:.2f}%,{self.percent_coverage:.2f}%'
#         s += ",partial" if self.partial else ""
#         s += ',truncated' if self.phenotype == "truncated" else ""
#         s += ",below_id_threshold" if self.below_threshold else ""
#         return s
#
#     @classmethod
#     def from_dict(cls, d: dict, **kwargs) -> GeneResult:
#         return cls(
#             id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'], dna_seq=Seq(d['dna_seq']),
#             protein_seq=Seq(d['protein_seq']), below_threshold=True if d['below_threshold'] == 'True' else False,
#             phenotype=d['phenotype'], gene_type=d['gene_type'], partial=True if d['partial'] == 'True' else False,
#             percent_identity=float(d['percent_identity']), percent_coverage=float(d['percent_coverage']), **kwargs
#         )
#
#     def format(self, format_spec) -> Union[str, dict]:
#         if format_spec == 'ffn':
#             if len(self.dna_seq) == 0:
#                 warn(f'No DNA sequence for {self}')
#                 return ""
#             return (f'>{self.gene.name} {self.result.genome}|{self.id}:{self.start}-{self.end}{self.strand}\n'
#                     f'{self.dna_seq}\n')
#         if format_spec == 'faa':
#             if len(self.protein_seq) == 0:
#                 warn(f'No protein sequence for {self.__repr__()}')
#                 return ""
#             return (f'>{self.gene.name} {self.result.genome}|{self.id}:{self.start}-{self.end}{self.strand}\n'
#                     f'{self.protein_seq}\n')
#         if format_spec == 'json':
#             return {
#                 'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
#                 'dna_seq': str(self.dna_seq), 'protein_seq': str(self.protein_seq), 'partial': str(self.partial),
#                 'below_threshold': str(self.below_threshold), 'phenotype': self.phenotype, 'gene_type': self.gene_type,
#                 'percent_identity': str(self.percent_identity), 'percent_coverage': str(self.percent_coverage),
#                 'gene': self.gene.name, 'piece': self.piece.__repr__() if self.piece else '',
#             }
#         raise ValueError(f"Unknown format specifier {format_spec}")
#
#
# class LocusFormatter:
#     def format(self, format_spec):
#         if format_spec == 'fna':
#             if len(self.seq) == 0:
#                 warn(f'No DNA sequence for {self}', UserWarning)
#                 return ""
#             return f'>{self.name}\n{self.seq}\n'
#         if format_spec in {'ffn', 'faa'}:
#             return ''.join([gene.format(format_spec) for gene in self])
#         raise ValueError(f'Invalid format specifier: {format_spec}')
#
#     def write(self, fna: Union[str, Path, IO[str]] = None, ffn: Union[str, Path, IO[str]] = None,
#               faa: Union[str, Path, IO[str]] = None):
#         """Write the serotyping result to files or file handles."""
#         for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
#             if f:
#                 if isinstance(f, TextIOBase):
#                     f.write(self.format(fmt))
#                 elif isinstance(f, (Path, str)):
#                     Path(f / f'{self.name.replace("/", "_")}.{fmt}').write_text(self.format(fmt))
#
#
#
# class TabularWriter(ABC):
#     __slots__ = ('_handle', '_handle_path', '_output_header')
#     _DELIMITER = '\t'
#     _HEADER = ''
#
#     def __init__(self, handle: IO[str], output_header: bool = True) -> None:
#         self._handle = handle
#         self._handle_path = None
#         if handle_name := getattr(self._handle, 'name', None):
#             if handle_name != '<stdout>':
#                 self._handle_path = Path(handle_name)
#         self._output_header = output_header
#
#     def _header_written(self) -> bool:
#         return (not self._output_header) and self._handle_path and self._handle_path.stat().st_size
#
#     def write_header(self):
#         """Write appropriate header to a file handle."""
#         self._handle.write(self._HEADER)
#         self._output_header = False
#
#     @abstractmethod
#     def write_one(self):
#         ...
#
#
# class ResultTabularWriter(TabularWriter):
#     _HEADER = ('Assembly\tBest match locus\tBest match type\tMatch confidence\tProblems\tIdentity\tCoverage\t'
#                'Length discrepancy\tExpected genes in locus\tExpected genes in locus, details\t'
#                'Missing expected genes\tOther genes in locus\tOther genes in locus, details\t'
#                'Expected genes outside locus\tExpected genes outside locus, details\t'
#                'Other genes outside locus\tOther genes outside locus, details\t'
#                'Truncated genes, details\tExtra genes, details\n')
#
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#
#     def _write_one(self, result: SerotyperResult) -> None:
#         out = self._DELIMITER.join(
#             [
#                 self.sample_name,
#                 self.best_match.name,
#                 self.phenotype,
#                 self.confidence,
#                 self.problems,
#                 f"{self.percent_identity:.2f}%",
#                 f"{self.percent_coverage:.2f}%",
#                 f"{self.__len__() - len(self.best_match)} bp" if len(self.pieces) == 1 else 'n/a',
#                 f"{(n_inside := len({i.gene.name for i in self.expected_genes_inside_locus}))} / {(n_expected := len(self.best_match.seqs))} ({100 * n_inside / n_expected:.2f}%)",
#                 ';'.join(map(str, self.expected_genes_inside_locus)),
#                 ';'.join(self.missing_genes),
#                 f"{len(self.unexpected_genes_inside_locus)}",
#                 ';'.join(map(str, self.unexpected_genes_inside_locus)),
#                 f"{(n_outside := len({i.gene.name for i in self.expected_genes_outside_locus}))} / {n_expected} ({100 * n_outside / n_expected:.2f}%)",
#                 ';'.join(map(str, self.expected_genes_outside_locus)),
#                 f"{len(self.unexpected_genes_outside_locus)}",
#                 ';'.join(map(str, self.unexpected_genes_outside_locus)),
#                 ';'.join(map(str, filter(lambda z: z.consequence == "truncated", self))),
#                 ';'.join(map(str, self.extra_genes))
#             ]
#         ) + "\n"
#
#
# class ResultSerializer:
#     _RESULT_ATTRS = []
#     _PEICE_ATTRS = []
#     _GENE_ATTRS = []
#
#     __slots__ = ('_handle',)
#     _DELIMITER = '\t'
#     _HEADER = ''
#
#     def __init__(self, handle: IO[str]):
#         self._handle = handle
#
#     def _write_one(self, result: SerotyperResult):
#         return json_dumps(
#         ) + "\n"
#
#     def _serialize_result(self, result: SerotyperResult) -> dict:
#         return {
#             'sample_name': result.genome, 'best_match': result.best_match_locus.name, 'confidence': result.confidence,
#             'phenotype': result.phenotype, 'problems': result.problems,
#             'percent_identity': str(result.percent_identity),
#             'percent_coverage': str(result.percent_coverage), 'missing_genes': result.missing_genes
#         } | {
#             attr: [i.format(format_spec) for i in getattr(result, attr)] for attr in {
#                 'pieces', 'expected_genes_inside_locus', 'unexpected_genes_inside_locus',
#                 'expected_genes_outside_locus', 'unexpected_genes_outside_locus', 'extra_genes'
#             }
#         }
#
#     def _serialize_piece(self, locus_piece: Record) -> dict: ...
#
#     def _serialize_gene(self, gene_result: GeneResult) -> dict:
#         return {
#             'id': gene_result.id, 'start': str(self.start), 'end': str(gene_result.end), 'strand': gene_result.strand,
#             'dna_seq': str(gene_result.dna_seq), 'protein_seq': str(gene_result.protein_seq),
#             'partial': str(gene_result.partial),
#             'below_threshold': str(gene_result.below_threshold), 'phenotype': gene_result.consequence,
#             'gene_type': gene_result.gene_type,
#             'percent_identity': str(gene_result.percent_identity),
#             'percent_coverage': str(gene_result.percent_coverage),
#             'gene': gene_result.gene.name, 'piece': gene_result.piece.__repr__() if gene_result.piece else '',
#         }
#
#
# class ScoresTabularWriter(TabularWriter):
#     _HEADER = 'Assembly\tLocus\tAS\tmlen\tblen\tq_len\tgenes_found\tgenes_expected\n'
#
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#
#
# class ResultDeserializer:
#     def __init__(self, handle: IO[str], database: Database):
#         self._handle = handle
#         self._database = database
#
#     def _read_one(self) -> 'SerotyperResult':
#         if not (best_match := db.loci.get(d['best_match'])):
#             raise TypingResultError(f"Best match {d['best_match']} not found in database")
#         result = SerotyperResult(sample_name=d['sample_name'], db=db, best_match=best_match,
#                                  missing_genes=d['missing_genes'])
#         # Set the cached properties
#         result._percent_identity = float(d['percent_identity'])
#         result._percent_coverage = float(d['percent_coverage'])
#         result._phenotype = d['phenotype']
#         result._problems = d['problems']
#         result._confidence = d['confidence']
#         # Add the pieces and create the gene results
#         result.pieces = [Record.from_dict(i, result=result) for i in d['pieces']]
#         pieces = {i.__repr__(): i for i in result.pieces}
#         gene_results = {}  # This was previously a dict comp, but we need to check the gene is in the database, see #31
#         for r in chain(d['expected_genes_inside_locus'], d['unexpected_genes_inside_locus'],
#                        d['expected_genes_outside_locus'], d['unexpected_genes_outside_locus'],
#                        d['extra_genes']):
#             if not (gene := db.seqs.get(r['gene'])) and not (gene := db.extra_genes.get(r['gene'])):
#                 raise TypingResultError(f"Gene {r['gene']} not found in database")
#             x = GeneResult.from_dict(r, result=result, piece=pieces.get(r['piece']), gene=gene)
#             gene_results[x.__repr__()] = x
#
#         for gene_result in gene_results.values():
#             self.add_gene_result(gene_result)
#         return self
#
#     def _read_gene(self, d: dict) -> 'GeneResult':
#         return GeneResult(
#             id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'], dna_seq=Seq(d['dna_seq']),
#             protein_seq=Seq(d['protein_seq']), below_threshold=True if d['below_threshold'] == 'True' else False,
#             phenotype=d['phenotype'], gene_type=d['gene_type'], partial=True if d['partial'] == 'True' else False,
#             percent_identity=float(d['percent_identity']), percent_coverage=float(d['percent_coverage']), **kwargs
#         )
