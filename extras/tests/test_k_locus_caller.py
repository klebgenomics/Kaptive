#!/usr/bin/env python
'''
K locus caller tests

This script executes some unit tests for the K locus caller script.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import unittest
import sys
import os


# Import the script to be tested, whether this test script is being run from the tests directory or
# the parent directory with the main script.
test_dir = os.path.dirname(os.path.realpath(__file__))
script_dir = os.getcwd()
while not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
    script_dir = os.path.dirname(script_dir)
if not os.path.isfile(os.path.join(script_dir, 'k_locus_caller.py')):
    print('Could not find k_locus_caller.py')
    sys.exit(1)
sys.path.append(script_dir)
from k_locus_caller import *

def delete_blast_databases():
    '''
    Delete any BLAST databases that might exist in the test directory
    '''
    for filename in os.listdir(test_dir):
        if filename.endswith('.nin') or filename.endswith('.nhr') or filename.endswith('.nsq'):
            os.remove(os.path.join(test_dir, filename))



class TestIntRange(unittest.TestCase):

    def test_constructor(self):
        int_range = IntRange([(1, 10), (11, 20), (15, 30)])
        self.assertEqual(int_range.ranges, [(1, 10), (11, 30)])
        int_range = IntRange([(1, 10), (20, 20)])
        self.assertEqual(int_range.ranges, [(1, 10)])

    def test_add_range(self):
        int_range = IntRange()
        self.assertEqual(int_range.ranges, [])
        int_range.add_range(5, 10)
        self.assertEqual(int_range.ranges, [(5, 10)])
        int_range.add_range(8, 20)
        self.assertEqual(int_range.ranges, [(5, 20)])
        int_range.add_range(12, 2)
        self.assertEqual(int_range.ranges, [(2, 20)])
        int_range.add_range(30, 40)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 40)])
        int_range.add_range(40, 50)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 50)])
        int_range.add_range(51, 60)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 50), (51, 60)])
        int_range.add_range(50, 50)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 50), (51, 60)])
        int_range.add_range(51, 51)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 50), (51, 60)])
        int_range.add_range(50, 51)
        self.assertEqual(int_range.ranges, [(2, 20), (30, 60)])
        int_range.add_range(20, 29)
        self.assertEqual(int_range.ranges, [(2, 29), (30, 60)])
        int_range.add_range(20, 30)
        self.assertEqual(int_range.ranges, [(2, 60)])

    def test_add_ranges(self):
        int_range = IntRange()
        self.assertEqual(int_range.ranges, [])
        int_range.add_ranges([(1, 10)])
        self.assertEqual(int_range.ranges, [(1, 10)])
        int_range.add_ranges([(40, 50), (5, 15)])
        self.assertEqual(int_range.ranges, [(1, 15), (40, 50)])
        int_range.add_ranges([(55, 45), (-3, -10)])
        self.assertEqual(int_range.ranges, [(-10, -3), (1, 15), (40, 55)])
        int_range.add_ranges([(-3, 1), (15, 40)])
        self.assertEqual(int_range.ranges, [(-10, 55)])

    def test_contains(self):
        int_range_1 = IntRange([(5, 10), (20, 30)])
        self.assertTrue(int_range_1.contains(IntRange([(5, 10), (20, 30)])))
        self.assertTrue(int_range_1.contains(IntRange([(6, 8), (9, 10), (5, 6), (22, 23), (26, 28)])))
        self.assertFalse(int_range_1.contains(IntRange([(6, 11)])))
        self.assertFalse(int_range_1.contains(IntRange([(4, 10)])))
        self.assertFalse(int_range_1.contains(IntRange([(1, 3)])))
        self.assertFalse(int_range_1.contains(IntRange([(40, 50)])))
        self.assertFalse(int_range_1.contains(IntRange([(8, 25)])))

    def test_get_total_length(self):
        self.assertEqual(IntRange([(0, 10)]).get_total_length(), 10)
        self.assertEqual(IntRange([(0, 10), (30, 40)]).get_total_length(), 20)
        self.assertEqual(IntRange([(0, 10), (5, 15), (10, 20)]).get_total_length(), 20)
        self.assertEqual(IntRange([(0, 10), (5, 15), (10, 20), (100, 100), (200, 200)]).get_total_length(), 20)



class TestAssembly(unittest.TestCase):

    def test_constructor(self):
        delete_blast_databases()
        assembly_spades = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        assembly_other = Assembly(os.path.join(test_dir, 'test_assembly_other_format.fasta'))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_spades_format.fasta.nin')))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_spades_format.fasta.nhr')))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_spades_format.fasta.nsq')))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_other_format.fasta.nin')))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_other_format.fasta.nhr')))
        self.assertTrue(os.path.isfile(os.path.join(test_dir, 'test_assembly_other_format.fasta.nsq')))
        self.assertEqual(assembly_spades.fasta, os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        self.assertEqual(assembly_other.fasta, os.path.join(test_dir, 'test_assembly_other_format.fasta'))
        self.assertEqual(len(assembly_spades.contigs), 4)
        self.assertEqual(assembly_spades.name, 'test_assembly_spades_format')
        self.assertEqual(assembly_other.name, 'test_assembly_other_format')
        self.assertEqual(len(assembly_spades.contigs), 4)
        self.assertEqual(len(assembly_other.contigs), 4)
        self.assertTrue('contig_1' in assembly_other.contigs)
        self.assertTrue('contig_2' in assembly_other.contigs)
        self.assertTrue('contig_2 extra header parts' not in assembly_other.contigs)
        delete_blast_databases()



class TestAssemblyPiece(unittest.TestCase):

    def test_constructor(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        self.assertEqual(piece.assembly.name, 'test_assembly_spades_format')
        self.assertEqual(piece.contig_name, 'NODE_5_length_150905_cov_4.42519_ID_15485')
        self.assertEqual(piece.start, 376)
        self.assertEqual(piece.end, 684)
        self.assertEqual(piece.strand, '+')
        self.assertEqual(piece.blast_hits, [])
        delete_blast_databases()

    def test_get_header(self):
        assembly_spades = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        assembly_other = Assembly(os.path.join(test_dir, 'test_assembly_other_format.fasta'))
        spades_piece = AssemblyPiece(assembly_spades, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        other_piece = AssemblyPiece(assembly_other, 'contig_1', 376, 684, '+')
        self.assertEqual(spades_piece.get_header(), 'NODE_5_377_to_684_+_strand')
        self.assertEqual(other_piece.get_header(), 'contig_1_377_to_684_+_strand')
        delete_blast_databases()

    def test_get_sequence(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        self.assertEqual(piece.get_sequence(), 'ACCGCCGTCGGTCATCAGGCGGATGCCGCCGGTAGCGCACAGCGGCCCTTCATCGACGGCGGCGACGGCGAACTGGTACCAGTTGTCGCCCATTTCCGGGCGCGGCAGGCCGAAGATCTCCACGGCTACGCCAAACTCAAAGGTGCACAGGCCATCGTAGGCGAGGGCGACCACCTGATGACGTGACAGTGTCACAGAAGGGGTTGTCAGGATGTTAGCGTTTTCTGTCATGGTTATCGGCAGATCAGGGTGATGGGCGTCAGTAAAGTAGTGTTATCACAGACAAGGAGAAGACGCCATGAGTGTAG')
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '-')
        self.assertEqual(piece.get_sequence(), 'CTACACTCATGGCGTCTTCTCCTTGTCTGTGATAACACTACTTTACTGACGCCCATCACCCTGATCTGCCGATAACCATGACAGAAAACGCTAACATCCTGACAACCCCTTCTGTGACACTGTCACGTCATCAGGTGGTCGCCCTCGCCTACGATGGCCTGTGCACCTTTGAGTTTGGCGTAGCCGTGGAGATCTTCGGCCTGCCGCGCCCGGAAATGGGCGACAACTGGTACCAGTTCGCCGTCGCCGCCGTCGATGAAGGGCCGCTGTGCGCTACCGGCGGCATCCGCCTGATGACCGACGGCGGT')
        delete_blast_databases()

    def test_get_bandage_range(self):
        assembly_spades = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        assembly_other = Assembly(os.path.join(test_dir, 'test_assembly_other_format.fasta'))
        spades_piece = AssemblyPiece(assembly_spades, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        other_piece = AssemblyPiece(assembly_other, 'contig_1', 376, 684, '+')
        self.assertEqual(spades_piece.get_bandage_range(), '(377) 5+ (684)')
        self.assertEqual(other_piece.get_bandage_range(), '(377) contig_1+ (684)')
        delete_blast_databases()

    def test_get_length(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        self.assertEqual(piece.get_length(), 308)
        delete_blast_databases()

    def test_get_sequence_short(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        self.assertEqual(piece.get_sequence_short(), 'ACCGCC...GTGTAG (308 bp)')
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 380, '+')
        self.assertEqual(piece.get_sequence_short(), 'ACCG (4 bp)')
        delete_blast_databases()

    def test_combine(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece_1 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 690, 700, '+')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece, None)
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 685, 700, '+')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece, None)
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 684, 700, '+')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece.get_length(), 324)
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 200, 700, '+')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece.get_length(), 500)
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 200, 700, '-')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece, None)
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 100, 500, '+')
        combined_piece = piece_1.combine(piece_2)
        self.assertEqual(combined_piece.get_length(), 584)
        self.assertEqual(combined_piece.get_sequence(), 'CAGCAAAGCATCGCCCACGTAAAGGACGTCCTCCACCACCTGGATCTGCGGAAAACGGGACTGCAGCGCCGCGGTGTAGCGCCAGTGGGTTGTCGCCTGGCGGCCATTGAGCAGCCCGGCGGCGGCGAGAACGAACACCCCGGAACAGATAGAAATAATCCGACATCCGCGGGCGTGGGCGGAAGCCAGCGCGGCGCACAGCGCCTCGGGGACCGGGGCGTCGACGCCGCGCCAGCCGGGGACCACGATCGTATCCGCCTGGGCGAGAAGTTCAGGACCGCCGTCGGTCATCAGGCGGATGCCGCCGGTAGCGCACAGCGGCCCTTCATCGACGGCGGCGACGGCGAACTGGTACCAGTTGTCGCCCATTTCCGGGCGCGGCAGGCCGAAGATCTCCACGGCTACGCCAAACTCAAAGGTGCACAGGCCATCGTAGGCGAGGGCGACCACCTGATGACGTGACAGTGTCACAGAAGGGGTTGTCAGGATGTTAGCGTTTTCTGTCATGGTTATCGGCAGATCAGGGTGATGGGCGTCAGTAAAGTAGTGTTATCACAGACAAGGAGAAGACGCCATGAGTGTAG')
        delete_blast_databases()

    def test_overlaps(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        self.assertTrue(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 100, 400))
        self.assertFalse(piece.overlaps('contig_1', 100, 400))
        self.assertTrue(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 100, 377))
        self.assertFalse(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 100, 376))
        self.assertTrue(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 400, 1000))
        self.assertTrue(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 683, 1000))
        self.assertFalse(piece.overlaps('NODE_5_length_150905_cov_4.42519_ID_15485', 684, 1000))
        delete_blast_databases()

    def test_hit_coordinates(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 376, 684, '+')
        piece.blast_hits.append(BlastHit('gene1\tNODE_5_length_150905_cov_4.42519_ID_15485\t10\t90\t400\t480\t2e-40\t12345\t81\t99.4\t81\tACGT'))
        piece.blast_hits.append(BlastHit('gene1\tNODE_5_length_150905_cov_4.42519_ID_15485\t5\t95\t400\t480\t2e-40\t12345\t81\t99.4\t81\tACGT'))
        piece.blast_hits.append(BlastHit('gene1\tNODE_5_length_150905_cov_4.42519_ID_15485\t200\t290\t400\t480\t2e-40\t12345\t81\t99.4\t81\tACGT'))
        piece.blast_hits.append(BlastHit('gene1\tNODE_5_length_150905_cov_4.42519_ID_15485\t100\t190\t400\t480\t2e-40\t12345\t81\t99.4\t81\tACGT'))
        self.assertEqual(piece.earliest_hit_coordinate(), 4)
        self.assertEqual(piece.latest_hit_coordinate(), 290)
        delete_blast_databases()



class TestKLocus(unittest.TestCase):

    def get_k1_seq(self):
        k_refs = load_fasta(os.path.join(test_dir, 'test_k_refs.fasta'))
        return k_refs[0][1]

    def test_constructor(self):
        k1 = KLocus('K1', self.get_k1_seq(), os.path.join(test_dir, 'test_gene_list.txt'))
        self.assertEqual(k1.name, 'K1')
        self.assertTrue(k1.seq.startswith('atgaatatggcgaatttgaaagcggtt'))
        self.assertEqual(len(k1.gene_names), 20)
        self.assertEqual(k1.gene_names[0], '36__K1__K1-CDS1-galF__00064')
        self.assertEqual(k1.gene_names[19], '18__K1__K1-CDS9-wzy__00060')

        

class TestFunctions(unittest.TestCase):

    def test_fill_assembly_piece_gaps(self):
        assembly = Assembly(os.path.join(test_dir, 'test_assembly_spades_format.fasta'))
        piece_1 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 300, 600, '+')
        piece_2 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 700, 900, '+')
        piece_3 = AssemblyPiece(assembly, 'NODE_5_length_150905_cov_4.42519_ID_15485', 950, 1000, '+')
        gap_filled_pieces = fill_assembly_piece_gaps([piece_1, piece_2, piece_3], 10)
        self.assertEqual(len(gap_filled_pieces), 3)
        gap_filled_pieces = fill_assembly_piece_gaps([piece_1, piece_2, piece_3], 49)
        self.assertEqual(len(gap_filled_pieces), 3)
        gap_filled_pieces = fill_assembly_piece_gaps([piece_1, piece_2, piece_3], 50)
        self.assertEqual(len(gap_filled_pieces), 2)
        self.assertEqual(gap_filled_pieces[0].start, 300)
        self.assertEqual(gap_filled_pieces[0].end, 600)
        self.assertEqual(gap_filled_pieces[1].start, 700)
        self.assertEqual(gap_filled_pieces[1].end, 1000)
        gap_filled_pieces = fill_assembly_piece_gaps([piece_1, piece_2, piece_3], 99)
        self.assertEqual(len(gap_filled_pieces), 2)
        gap_filled_pieces = fill_assembly_piece_gaps([piece_1, piece_2, piece_3], 100)
        self.assertEqual(len(gap_filled_pieces), 1)
        self.assertEqual(gap_filled_pieces[0].start, 300)
        self.assertEqual(gap_filled_pieces[0].end, 1000)
        delete_blast_databases()



if __name__ == '__main__':
    unittest.main()