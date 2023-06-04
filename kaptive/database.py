"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/kaptive-mapper

This file is part of kaptive-mapper. kaptive-mapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. kaptive-mapper is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with kaptive-mapper.
If not, see <https://www.gnu.org/licenses/>.
"""
from pathlib import Path
from Bio import SeqIO
from collections import OrderedDict

from .log import quit_with_error
from .misc import quick_translate


class Database(object):
    def __init__(self, kaptive_refs: Path, locus_label: str, type_label: str):
        self.kaptive_refs = kaptive_refs
        self.locus_label = locus_label
        self.type_label = type_label
        self.loci = self.parse_genbank()
        self.ref_types = {
            locus.name: locus.type for locus in self.loci.values() if not locus.name.startswith('Extra_genes')
        }
        self.genes = {gene.name: gene for locus in self.loci.values() for gene in locus.genes}
        self.special_logic = self.load_special_logic() if any(
            t == 'special logic' for t in self.ref_types.values()) else []

    def __repr__(self):
        return self.kaptive_refs.stem

    def __str__(self):
        return self.kaptive_refs.stem

    def __len__(self):
        return len(self.loci)

    def get_fasta(self):
        return ''.join(locus.get_fasta() for locus in self.loci.values())

    def get_gene_fasta(self, translation=False):
        return ''.join(gene.get_fasta(translation) for gene in self.genes.values())

    def get_gff(self):
        return ''.join(locus.get_gff() for locus in self.loci.values())

    def get_bed(self):
        return ''.join(locus.get_bed() for locus in self.loci.values())

    def parse_genbank(self):
        refs = {}
        if self.kaptive_refs:
            locus_label = check_label(
                self.kaptive_refs, self.locus_label
            ) if self.locus_label else find_label(self.kaptive_refs, 'locus')
            type_label = check_label(
                self.kaptive_refs, self.type_label
            ) if self.type_label else find_label(self.kaptive_refs, 'type', required=False)
            for record in SeqIO.parse(self.kaptive_refs, 'genbank'):
                locus_name, type_name = '', ''
                for feature in record.features:
                    if feature.type == 'source' and 'note' in feature.qualifiers:
                        for note in feature.qualifiers['note']:
                            if note.startswith(locus_label):
                                locus_name = get_name_from_note(note, locus_label)
                            elif note.startswith('Extra genes'):
                                locus_name = note.replace(':', '').replace(' ', '_')
                            elif type_label is not None and note.startswith(type_label):
                                type_name = get_name_from_note(note, type_label)

                if locus_name in refs:
                    quit_with_error('Duplicate reference locus name: ' + locus_name)

                gene_num, genes = 1, []
                for feature in record.features:
                    if feature.type == 'CDS':
                        genes.append(Gene(locus_name, gene_num, feature, record.seq))
                        gene_num += 1
                refs[locus_name] = Locus(locus_name, type_name, str(record.seq), genes, self)
                for gene in genes:  # Add the locus object to its respective genes, makes things easier downstream
                    gene.locus = refs[locus_name]
        return refs

    def load_special_logic(self) -> list:
        """
        This function loads special logic file if  any of the reference loci have a type of 'special logic',
        that implies that a corresponding file exists to describe that logic.
        Logic files are tab-delimited files with 3 columns: Locus name, Extra genes, Type.
        Returns a dictionary with the locus name as the key and the extra genes and type as the value.
        """
        special_logic_filename = self.kaptive_refs.with_suffix(".logic")
        if not special_logic_filename.exists():
            quit_with_error(f'No special logic file found for {self.kaptive_refs}')

        with open(special_logic_filename, 'rt') as f:
            lines = [[x.strip() for x in i.split('\t')] for i in f.read().splitlines()]
        if len(lines[0]) != 3:
            quit_with_error(f'Invalid special logic file found for {self.kaptive_refs}')

        logic = []
        for line in lines[1:]:
            logic_dict = {
                'locus': line[0], 'extra_loci': line[1].split(',') if line[1].lower() != 'none' else [], 'type': line[2]
            }
            if line[0] not in self.loci:
                quit_with_error(
                    f'Invalid locus name found in special logic file for {self.kaptive_refs}: {line[0]}')
            else:
                self.loci[line[0]].special_logic.append(logic_dict)
        return logic


class Locus(object):
    def __init__(self, name, type_name, seq, genes, db: 'Database'):
        self.name = name
        self.db = db
        self.type = type_name
        self.seq = seq
        self.genes = genes
        self.gene_names = [x.name for x in genes]
        self.special_logic = []

    def __repr__(self):
        return f'Locus {self.name}'

    def __len__(self):
        """Returns the locus sequence length."""
        return len(self.seq)

    def get_fasta(self):
        """Returns the locus sequence as a FASTA string."""
        return f'>{self.name}\n{self.seq}\n'

    def get_gene_fasta(self, translation=False):
        return ''.join(gene.get_fasta(translation) for gene in self.genes)

    def get_gff(self):
        return ''.join(gene.get_gff() for gene in self.genes)

    def get_bed(self):
        return ''.join(gene.get_bed() for gene in self.genes)


class Gene(object):
    """This class prepares and stores a gene taken from the input Genbank file."""

    def __init__(self, locus_name, num, feature, k_locus_seq):
        self.locus_name = locus_name
        self.feature = feature
        self.position_in_locus = num
        gene_num_string = str(num).zfill(2)
        self.name = locus_name + '_' + gene_num_string
        if 'gene' in feature.qualifiers:
            self.gene_name = feature.qualifiers['gene'][0]
            self.name += '_' + self.gene_name
        else:
            self.gene_name = None
        if 'product' in feature.qualifiers:
            self.product = feature.qualifiers['product'][0]
        else:
            self.product = None
        if 'EC_number' in feature.qualifiers:
            self.ec_number = feature.qualifiers['EC_number'][0]
        else:
            self.ec_number = None
        self.nuc_seq = feature.extract(k_locus_seq)
        self.prot_seq = quick_translate(self.nuc_seq, to_stop=True)
        self.nuc_seq = str(self.nuc_seq)
        self.locus = None

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.nuc_seq)

    def get_fasta(self, translation=False):
        """Returns the locus protein or nucleotide sequence as a FASTA string."""
        return f'>{self.name}\n{self.prot_seq if translation else self.nuc_seq}\n'

    def get_gff(self):
        """
        Returns the gene GFF strings compatible with bcftools csq.
        See format specification here: https://samtools.github.io/bcftools/bcftools.html#csq
        """

        strand = '+' if self.feature.strand == 1 else '-'

        gene_line = f"{self.locus_name}\tKaptive\tgene\t{self.feature.location.start + 1}\t" \
                    f"{self.feature.location.end}\t.\t{strand}\t.\t" \
                    f"ID=gene:{self.name};biotype=protein_coding;Name={self.gene_name};description={self.product};"

        mrna_line = f"{self.locus_name}\tKaptive\ttranscript\t{self.feature.location.start + 1}\t" \
                    f"{self.feature.location.end}\t.\t{strand}\t.\t" \
                    f"ID=transcript:{self.name}.t1;Parent=gene:{self.name};biotype=protein_coding;"

        cds_line = f"{self.locus_name}\tKaptive\tCDS\t{self.feature.location.start + 1}\t" \
                   f"{self.feature.location.end}\t.\t{strand}\t0\t" \
                   f"ID=CDS:{self.name}.cds1;Parent=transcript:{self.name}.t1;"

        return f"{gene_line}\n{mrna_line}\n{cds_line}\n"

    def get_bed(self):
        """Returns the gene BED string compatible with bcftools csq."""
        strand = '+' if self.feature.strand == 1 else '-'
        return f"{self.locus_name}\t{self.feature.location.start}\t{self.feature.location.end}\t" \
               f"{self.name}\t0\t{strand}\n"

    def get_reference_info_json_dict(self):
        reference_dict = OrderedDict()
        if self.gene_name:
            reference_dict['Gene'] = self.gene_name
        if self.product:
            reference_dict['Product'] = self.product
        if self.ec_number:
            reference_dict['EC number'] = self.ec_number
        reference_dict['Nucleotide length'] = len(self.nuc_seq)
        reference_dict['Protein length'] = len(self.prot_seq)
        reference_dict['Nucleotide sequence'] = self.nuc_seq
        reference_dict['Protein sequence'] = self.prot_seq
        return reference_dict


def check_label(genbank: Path, label: str):
    """
    Makes sure that every record in the Genbank file contains a note in the source feature
    beginning with the given label.
    """
    for record in SeqIO.parse(genbank, 'genbank'):
        found_label = False
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if note.startswith(label):
                        locus_name = get_name_from_note(note, label)
                        if locus_name:
                            found_label = True
        if not found_label:
            quit_with_error(f'{record.name} is missing a label\n'
                            f'The source feature must have a note qualifier '
                            f'beginning with "{label}:" followed by the relevant info')


def get_name_from_note(full_note, locus_label):
    """
    Extracts the part of the note following the label (and any colons, spaces or equals signs).
    """
    locus_name = full_note[len(locus_label):].strip()
    while locus_name.startswith(':') or locus_name.startswith(' ') or \
            locus_name.startswith('='):
        locus_name = locus_name[1:]
    return locus_name


def find_label(genbank: Path, text: str, required=True):
    """
    Automatically finds the label in the Genbank file which contains the specified text. For
    example, if the text is 'locus', then the Genbank file must have exactly one possible label
    containing 'locus' that is present in a note qualifier in the source feature for every record.
    If not, Kaptive will quit with an error.
    """
    possible_locus_labels = set()
    for record in SeqIO.parse(genbank, 'genbank'):
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if ':' in note:
                        note = note.split(':')[0].strip()
                        if text in note:
                            possible_locus_labels.add(note)
    if not possible_locus_labels:
        if required:
            quit_with_error(f'None of the records contain a valid {text} label')
        else:
            return None
    available_locus_labels = possible_locus_labels.copy()
    for record in SeqIO.parse(genbank, 'genbank'):
        locus_labels = set()
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if ':' in note:
                        locus_labels.add(note.split(':')[0].strip())
        if any(x == 'Extra genes' for x in locus_labels):
            continue
        if not locus_labels:
            quit_with_error(f'no possible {text} labels were found for {record.name}')
        previous_labels = available_locus_labels.copy()
        available_locus_labels = available_locus_labels.intersection(locus_labels)
        if not available_locus_labels:
            quit_with_error(f'{record.name} does not have a {text} label matching the previous records\n'
                            f'Previous record labels: {", ".join(list(previous_labels))}\n'
                            f'Labels in {record.name}: {", ".join(list(locus_labels))}')
    if len(available_locus_labels) > 1:
        quit_with_error(f'multiple possible {text} labels were found: {", ".join(list(available_locus_labels))}\n'
                        f'Please use the --{text}_label option to specify which to use')
    return list(available_locus_labels)[0]
