#!/usr/bin/env python

import logging

LOG = logging.getLogger(__name__)
# LOG.setLevel(logging.INFO)

# Set logging defaults
logging.basicConfig(format=('%(levelname)s '
                            '%(filename)s '
                            '%(funcName)s() '
                            '%(message)s'), level=logging.DEBUG)

import sys
import time
import unittest

import utils

class TestApp(unittest.TestCase):

    def __init__(self, testCaseNames):
        ''' Each test case specified below gets it instance of a unit
            test object with independent self definitions.
        '''
        super(TestApp, self).__init__(testCaseNames)

        utils.DEBUG = True

    # End of TestApp.__init__()

    def setUp(self):
        ''' TC setup'''
        LOG.debug('=======================================================')
        LOG.debug('Start of new Test Case:')


    # End of TestApp.setUp()

    def tearDown(self):
        '''TC tear down'''
        LOG.debug('Completion of this Test Case')
        LOG.debug('=======================================================')

    # End of TestApp.tearDown()

    def test_KmersFromText(self):

        text = 'CCGGCGTTAG'
        k = 4

        results = utils.KmersFromText(text, k)

        print ("Results are:", results)

        expected_results = [
                'CCGG',
                'CGGC',
                'GGCG',
                'GCGT',
                'CGTT',
                'GTTA',
                'TTAG'
            ]

        self.assertEqual(results, expected_results)

    # End of test_KmersFromText()

    def test_RandomWeighted(self):

        probabilities = [0.2, 0.4, 0.2]
        result_counts = [0, 0, 0]
        loop_count = 5000
        expected_result_counts = [loop_count/4, loop_count/2, loop_count/4]
        for i in range(loop_count):
            result_counts[utils.RandomWeighted(probabilities)] += 1

        LOG.info("Result counts:")
        for i in range(len(result_counts)):
            print ("   Result for", i, "is", result_counts[i])

        for i in range(len(result_counts)):
            # Ensure results are within 10% of expected
            self.assertLess(abs(expected_result_counts[i] - result_counts[i]),
                            loop_count/10,
                            msg="Results out of range for instance " + str(i) +
                            " Actual result: " + str(result_counts[i]) +
                            " Expected result: " + str(expected_result_counts[i]))

    # End of test_RandomWeighted()

    def test_BaseToNumber(self):

        self.assertEqual(utils.BaseToNumber('A'), 0)
        self.assertEqual(utils.BaseToNumber('C'), 1)
        self.assertEqual(utils.BaseToNumber('G'), 2)
        self.assertEqual(utils.BaseToNumber('T'), 3)

    # End of test_BaseToNumber()

    def compare_results(self, true_results, expected_results):
        self.assertEqual(len(true_results), len(expected_results), msg =
                         "Unexpected lenght of results: Expected " +
                         str(len(expected_results)) + " Found " +
                         str(len(true_results)))

        for item in expected_results:
            self.assertTrue(item in true_results, msg =
                            str(item) + " expected, but not found in " +
                            str(true_results))
        for item in true_results:
            self.assertTrue(item in expected_results, msg =
                            str(item) + " found, but not expected in " +
                            str(true_results))

    # End of TestApp.compare_results()

    def test_MotifEnumeration(self):

        LOG.info("Test with test input")
        result = utils.MotifEnumeration2(
            3, 1, ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'])
        self.compare_results(result, ['ATA', 'ATT', 'GTT', 'TTT'])

        LOG.info("Test with params for stepik test #5")
        result = utils.MotifEnumeration2(
            3, 0, ['AAAAA', 'AAAAA', 'AACAA'])
        import pdb; pdb.set_trace()
        self.compare_results(result, [])

        LOG.info("Test with d=0")
        result = utils.MotifEnumeration2(
            3, 0, ['ATCTGGC', 'TATCTTA', 'CGGTATC', 'GAAATCT'])
        self.compare_results(result, ['ATC'])

        LOG.info("Test with k=0, d=0")
        result = utils.MotifEnumeration2(
            0, 0, ['ATCTGGC', 'TATCTTA', 'CGGTATC', 'GAAATCT'])
        self.compare_results(result, [''])

        LOG.info("Test with k=1, d=0")
        result = utils.MotifEnumeration2(
            1, 0, ['ATCTGGC', 'TATCTTA', 'CGGTATC', 'GAAATCT'])
        self.compare_results(result, ['A', 'C', 'T'])

        LOG.info("Test with k=1, d=1")
        result = utils.MotifEnumeration2(
            1, 1, ['ATCTGGC', 'TATCTTA', 'CGGTATC', 'GAAATCT'])
        self.compare_results(result, ['A', 'C', 'G', 'T'])

        LOG.info("Test with small strands, k = d = length of strand")
        result = utils.MotifEnumeration2(
            2, 2, ['AT', 'TA', 'CG', 'GA'])
        self.compare_results(result, ['AA', 'AC', 'AG', 'AT',
                                  'CA', 'CC', 'CG', 'CT',
                                  'GA', 'GC', 'GG', 'GT',
                                  'TA', 'TC', 'TG', 'TT'])

        LOG.info("Test with d = k-1")
        result = utils.MotifEnumeration2(
            3, 2, ['ATCTGGC', 'TATCTTA', 'CGGTATC', 'GAAATCT'])
        self.compare_results(result,
            ['AAA', 'AAC', 'AAG', 'AAT',
             'ACA', 'ACC', 'ACG', 'ACT',
             'AGA', 'AGC', 'AGG', 'AGT',
             'ATA', 'ATC', 'ATG', 'ATT',
             'CAA', 'CAC', 'CAG', 'CAT',
             'CCA', 'CCC', 'CCG', 'CCT',
             'CGA', 'CGC',        'CGT',
             'CTA', 'CTC', 'CTG', 'CTT',
             'GAA', 'GAC', 'GAG', 'GAT',
             'GCA', 'GCC', 'GCG', 'GCT',
             'GGA', 'GGC',        'GGT',
             'GTA', 'GTC', 'GTG', 'GTT',
             'TAA', 'TAC', 'TAG', 'TAT',
             'TCA', 'TCC', 'TCG', 'TCT',
             'TGA', 'TGC', 'TGG', 'TGT',
             'TTA', 'TTC', 'TTG', 'TTT'])

    # End of test_MotifEnumeration()

    def test_GenProfileFromMotifs(self):

        motifs = [
                ['T', 'A', 'A', 'C'],
                ['G', 'T', 'C', 'T'],
                ['A', 'C', 'T', 'A'],
                ['A', 'G', 'G', 'T']
            ]

        profile = utils.GenProfileFromMotifs(motifs)

        print ("Profile returned is:")
        for row in profile:
            print ("   ", row)

        expected_profile = [
                [3.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [1.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 1.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 3.0/8.0],
            ]

        self.assertEqual(profile, expected_profile)

    # End of test_GenProfileFromMotifs()

    def test_GenProfileProbabilitiesKmer(self):

        kmers = ['CCGG', 'CGGC', 'GGCG', 'GCGT', 'CGTT', 'GTTA', 'TTAG']
        profile = [
                [3.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [1.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 1.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 3.0/8.0],
            ]
        probabilities = []

        for kmer in kmers:
            probabilities.append(utils.GenProfileProbabilitiesKmer(kmer, profile))

        print ("Found:", probabilities)

        denom = 8.0**4.0
        expected_probabilities = [
                4/denom, 8/denom, 8/denom, 24/denom,
                12/denom, 16/denom, 8/denom
            ]

        self.assertEqual(probabilities, expected_probabilities)

    # End of test_GenProfileProbabilitiesKmer()

    def test_ProfileRandomlyGeneratedKmer(self):

        text = 'CCGGCGTTAG'
        k = 4
        profile = [
                [3.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [1.0/8.0, 2.0/8.0, 2.0/8.0, 2.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 1.0/8.0],
                [2.0/8.0, 2.0/8.0, 2.0/8.0, 3.0/8.0],
            ]

        results = {
                'CCGG' : 0,
                'CGGC' : 0,
                'GGCG' : 0,
                'GCGT' : 0,
                'CGTT' : 0,
                'GTTA' : 0,
                'TTAG' : 0
            }

        iter_count = 10000

        for j in range(iter_count):
            kmer = utils.ProfileRandomlyGeneratedKmer(text, k, profile)
            results[kmer] += 1

        # Check the results

        denom = 80.0
        expected_results = {
                'CCGG' : ( 4/denom) * iter_count,
                'CGGC' : ( 8/denom) * iter_count,
                'GGCG' : ( 8/denom) * iter_count,
                'GCGT' : (24/denom) * iter_count,
                'CGTT' : (12/denom) * iter_count,
                'GTTA' : (16/denom) * iter_count,
                'TTAG' : ( 8/denom) * iter_count
            }

        print ("  Kmer:   Results:  Expected Results:")
        for kmer in utils.KmersFromText(text, k):
            print ("   {kmer}   {result:6}       {expected_result:6}".format(
                  kmer=kmer, result=results[kmer],
                  expected_result=int(expected_results[kmer])))

        for kmer in sorted(expected_results.keys()):
            self.assertLess(abs(results[kmer] - int(expected_results[kmer])),
                            iter_count/50,
                            msg="Results out of range for '" + kmer +
                            "'. Actual result: " + str(results[kmer]) +
                            ", Expected result: " + str(expected_results[kmer]))

    # End of test_ProfileRandomlyGeneratedKmer()

    def test_FindMostFrequentUsingDict(self):

        print ("test_FindMostFrequentUsingDict: Test edge conditions")
        result = utils.FindMostFrequentUsingDict("ABCDABCDAB", 2, 3)
        self.assertEqual(result, ['AB'])

        print ("test_FindMostFrequentUsingDict: Test overlap and no 't'")
        result = utils.FindMostFrequentUsingDict("ABCDABCD", 2)
        self.assertEqual(result, ['AB', 'BC', 'CD'])

        print ("test_FindMostFrequentUsingDict: Test no results")
        result = utils.FindMostFrequentUsingDict("ABCDABCD", 2, 3)
        self.assertEqual(result, [])

        print ("test_FindMostFrequentUsingDict: Test differing counts")
        result = utils.FindMostFrequentUsingDict("ABCDABABCDCDAB", 2, 3)
        self.assertEqual(result, ['AB', 'CD'])

        print ("test_FindMostFrequentUsingDict: Test differing counts, no 't'")
        result = utils.FindMostFrequentUsingDict("ABCDABABCDCDAB", 2)
        self.assertEqual(result, ['AB'])

    # End of test_FindMostFrequentUsingDict()

    def test_FindPatternClump_debug_dataset(self):
        print ("test_FindPatternClump: Debug dataset")
        result = utils.FindPatternClump(
            *get_FindPatternClump_data_from_file("data/ba1e_debug_dataset.txt"))
        self.assertEqual(result, ['CGACA', 'GAAGA'])
    # End of test_FindPatternClump_debug_dataset()

    def test_FindPatternClump_test_dataset(self):
        print ("test_FindPatternClump: Test dataset")
        result = utils.FindPatternClump(
            *get_FindPatternClump_data_from_file("data/ba1e_test_dataset.txt"))
        self.assertEqual(result, ['AAACCAGGTGG'])
    # End of test_FindPatternClump_test_dataset()

    def test_FindPatternClump_rosalind_final_dataset(self):
        print ("test_FindPatternClump: Test dataset")
        result = utils.FindPatternClump(
            *get_FindPatternClump_data_from_file("data/rosalind_ba1e.txt"))
        self.assertEqual(result, ['CCAACCCGTGTC', 'CGGACGGTACCC', 'GCTGACGGTAGA',
                                  'GTAGCGATTTAG', 'AGATTGAGTACA', 'CCTGAGACCGGT',
                                  'CTGAGACCGGTC', 'TGAGACCGGTCC', 'GAGACCGGTCCT',
                                  'AGACCGGTCCTG', 'GACCGGTCCTGA', 'ACCGGTCCTGAG',
                                  'CCGGTCCTGAGA', 'CGGTCCTGAGAC', 'GGTCCTGAGACC',
                                  'GTCCTGAGACCG', 'TCCTGAGACCGG',])
    # End of test_FindPatternClump_rosalind_final_dataset()

    def test_FindPatternClump_off_by_one_dataset(self):
        print ("test_FindPatternClump: Off-by-one dataset")
        result = utils.FindPatternClump(
            *get_FindPatternClump_data_from_file("data/ba1e_test_for_off_by_one.txt"))
        self.assertEqual(result, ['CGACA', 'GAAGA', 'TGTAA'])
    # End of test_FindPatternClump_off_by_one_dataset()

    def test_FindPatternClump_class_example(self):
        print ("test_FindPatternClump: class example")
        genome = "gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac".upper()
        result = utils.FindPatternClump(k=4, L=25, t=3, genome=genome)
        self.assertEqual(result, ['TGCA'])
    # End of test_FindPatternClump_class_example()

    def test_FindPatternClump_window_size(self):
        print ("test_FindPatternClump: Test window size")
        result = utils.FindPatternClump(2, 5, 2, "GCGATAATCGC")
        self.assertEqual(result, ['AT'])
        result = utils.FindPatternClump(2, 4, 2, "GCGATAATCGC")
        self.assertEqual(result, [])
        result = utils.FindPatternClump(2, 5, 2, "ATAATCGC")
        self.assertEqual(result, ['AT'])
        result = utils.FindPatternClump(2, 5, 2, "GCGATAAT")
        self.assertEqual(result, ['AT'])
        result = utils.FindPatternClump(2, 5, 2, "ATAAT")
        self.assertEqual(result, ['AT'])
    # End of test_FindPatternClump_window_size()

    def test_DeBruijnGraphs(self):

        print("\n<===== Testing DeBruijn Eulerian Cycle =====>\n")
        graph = utils.DeBruijnk(["0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2",
                                 "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6"])
        result = graph.EulerianCycle()
        string_result = '->'.join([str(x) for x in result])
        self.assertEqual(string_result, "6->5->4->2->1->0->3->2->6->8->7->9->6")

        print("\n<===== Testing DeBruijn Eulerian Path =====>\n")
        graph = utils.DeBruijnk(["0 -> 2", "1 -> 3", "2 -> 1", "3 -> 0,4",
                                 "6 -> 3,7", "7 -> 8", "8 -> 9", "9 -> 6"])
        result = graph.EulerianPath()
        string_result = '->'.join([str(x) for x in result])
        self.assertEqual(string_result, "6->7->8->9->6->3->0->2->1->3->4")

        print("\n<===== Testing DeBruijn Reconstruct Path String =====>\n")
        graph = utils.DeBruijnk(["CTTA", "ACCA", "TACC", "GGCT", "GCTT", "TTAC"])
        result = graph.ReconstructPathString()
        self.assertEqual(result, "GGCTTACCA")

        print("\n<===== Testing k-Universal string =====>\n")
        result = utils.KUniversalCircularString(4)
        self.assertEqual(result, '1110110010100001')

        print("\n<===== Testing MaximalNonBranchingPaths =====>\n")
        adjacency_list = ['1 -> 2', '2 -> 3', '3 -> 4,5', '6 -> 7', '7 -> 6']
        graph = utils.DeBruijnk(adjacency_list)
        result = graph.FindMaximalNonBranchingPaths()
        # result is a list of a list of edges, but we need to print out a list of nodes.
        node_paths = []
        for edge_path in result:
            node_path = [x.prev_node for x in edge_path]
            node_path.append(edge_path[-1].next_node)
            node_paths.append(node_path)
        str_result = [' -> '.join([x.label for x in node_path]) for node_path in node_paths]
        print("Result is:")
        for s_res in str_result:
            print("  ", s_res)
        self.assertTrue('1 -> 2 -> 3' in str_result, msg="'1 -> 2 -> 3' not in result")
        self.assertTrue('3 -> 4' in str_result, msg="'3 -> 4' not in result")
        self.assertTrue('3 -> 5' in str_result, msg="'3 -> 5' not in result")
        self.assertTrue('6 -> 7 -> 6' in str_result, msg="'6 -> 7 -> 6' not in result")

        print("\n<===== Testing Generating Contigs =====>\n")
        kmers = ['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']
        expected = ['AGA', 'ATG', 'ATG', 'CAT', 'GAT', 'TGGA', 'TGT']
        graph = utils.DeBruijnk(kmers)
        result = graph.GenContigs()
        print("Result is:", result)
        self.assertEqual(len(result), len(expected), msg="Unexpected length of results")
        for contig in result:
            self.assertTrue(contig in expected)

        print("\n<===== Testing StringSpelledByGappedPatters =====>\n")
        k = 4
        d = 2
        kdmers = ['GACC|GCGC', 'ACCG|CGCC', 'CCGA|GCCG', 'CGAG|CCGG', 'GAGC|CGGA']
        graph = utils.PairedDeBruijn(k,d,kdmers)
        result = graph.GetStringSpelledByGappedPatterns()
        print("Result is:", result)
        self.assertTrue(result == 'GACCGAGCGCCGGA', msg="Error, expected 'GACCGAGCGCCGGA'")

        print("\n<===== Testing AllGraphs Container =====>\n")
        # Input from coursera sample data
#        start_graph = utils.PairedDeBruijn(4, 2, [
#            "GAGA|TTGA", "TCGT|GATG", "CGTG|ATGT", "TGGT|TGAG", "GTGA|TGTT",
#            "GTGG|GTGA", "TGAG|GTTG", "GGTC|GAGA", "GTCG|AGAT"])
        # Input from book
        start_graph = utils.PairedDeBruijn(3, 1, [
            "TAA|GCC", "AAT|CCA", "CCA|GGG", "ATG|GAT", "GCC|TGG", "GGG|TGT",
            "ATG|CAT", "GGA|GTT", "CAT|GGA", "TGC|ATG", "TGG|ATG"])
        graph_container = utils.AllGraphs(start_graph)
        graph_container.SplitIntoSimpleDirectedGraphs()
        self.assertEqual(len(graph_container.all_graphs), 4)
        for graph in graph_container.all_graphs:
            path_string = graph.GetStringSpelledByGappedPatterns()
            print("Path string: ", path_string)

        print("\n<===== Testing Paired DeBruijn Reconstruct Path String =====>\n")
        graph = utils.PairedDeBruijn(4, 2, [
            "GAGA|TTGA", "TCGT|GATG", "CGTG|ATGT", "TGGT|TGAG", "GTGA|TGTT",
            "GTGG|GTGA", "TGAG|GTTG", "GGTC|GAGA", "GTCG|AGAT"])
        result = graph.ReconstructPathString()
        self.assertEqual(result, "GTGGTCGTGAGATGTTGA")

    # End of test_DeBruijnGraphs()

    def test_TranslateRNAToProtein(self):

        print("\n<===== Testing Translating RNA to Protein =====>\n")
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        protein = utils.TranslateRNAToProtein(rna)
        self.assertEqual(protein, 'MAMAPRTEINSTRING')

    # End of TestApp.test_TranslateRNAToProtein()

    def test_FindGeneForPeptide(self):

        print("\n<===== Testing Finding Gene for Peptide, using coursera sample data =====>\n")
        text = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        peptide = 'MA'

        results = utils.FindGeneForPeptide(text, peptide)

        expected_results = ['ATGGCC', 'GGCCAT', 'ATGGCC']

        self.assertEqual(len(results), len(expected_results))
        for res in results:
            self.assertTrue(res in expected_results)
        for res in expected_results:
            self.assertTrue(res in results)

        print("\n<===== Testing Finding Tyrocidin B1 gene in Bacillus brevis genome =====>\n")
        with open('Bacillus_brevis.txt', 'r') as f:
            temp = [x.strip() for x in f.readlines()]
            Bacillus_brevis_dna = ''.join(temp)
        Tyrocidine_B1 = 'VKLFPWFNQY'

        results = utils.FindGeneForPeptide(Bacillus_brevis_dna, Tyrocidine_B1)

        self.assertEqual(len(results), 0)

        print(results)

    # End of TestApp.test_FindGeneForPeptide()

    def test_GenTheoreticalSpectrum(self):

        print("\n<===== Testing GenTheoreticalSpectrum with cyclic coursera test input =====>\n")
        peptide = 'LEQN'
        result = utils.GenTheoreticalSpectrum(peptide, cyclic=True)
        expected_result = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        self.assertEqual(result, expected_result)

        print("\n<===== Testing GenTheoreticalSpectrum with linear coursera test input =====>\n")
        peptide = 'LEQN'
        result = utils.GenTheoreticalSpectrum(peptide, cyclic=False)
        expected_result = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
        self.assertEqual(result, expected_result)

    # End of TestApp.test_GenTheoreticalSpectrum()

    def test_CountingPeptidesWithGivenMass(self):

        print("\n<===== Testing CountingPeptidesWithGivenMass with single expected result =====>\n")
        test_mass = 57
        expected_result = 1
        result = utils.CountingPeptidesWithGivenMass(test_mass)
        self.assertEqual(result, expected_result)

        print("\n<===== Testing CountingPeptidesWithGivenMass with 3 matches =====>\n")
        test_mass = 128
        expected_result = 3
        result = utils.CountingPeptidesWithGivenMass(test_mass)
        self.assertEqual(result, expected_result)

        print("\n<===== Testing CountingPeptidesWithGivenMass with 17 matches =====>\n")
        test_mass = 241
        expected_result = 17
        import pdb; pdb.set_trace()
        result = utils.CountingPeptidesWithGivenMass(test_mass)
        self.assertEqual(result, expected_result)

        print("\n<===== Testing CountingPeptidesWithGivenMass with coursera test input =====>\n")
        import pdb; pdb.set_trace()
        test_mass = 1024
        expected_result = 14712706211
        result = utils.CountingPeptidesWithGivenMass(test_mass)
        self.assertEqual(result, expected_result)

    # End of test_CountingPeptidesWithGivenMass()

    def test_CyclopeptideSequencing(self):

        print("\n<===== Testing CyclopeptideSequencing with coursera test input =====>\n")
        spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
        expected_results = [
                [186, 128, 113],
                [186, 113, 128],
                [128, 186, 113],
                [128, 113, 186],
                [113, 186, 128],
                [113, 128, 186],
            ]
        results = utils.CyclopeptideSequencing(spectrum)
        self.assertEqual(len(results), len(expected_results))
        for peptide in results:
            self.assertTrue(peptide in expected_results)
            results.remove(peptide)
            expected_results.remove(peptide)

    # End of test_CyclopeptideSequencing()

# End of class TestApp

def get_FindPatternClump_data_from_file(fname):
    with open(fname, 'r') as f:
        genome = f.readline().strip()
        k, L, t = [int(x) for x in f.readline().split()]
    return k, L, t, genome


if __name__ == "__main__":

    LOG.debug("Unittest args: " + ' '.join(sys.argv))
    LOG.debug("Begining of the test suite:")
    LOG.debug("===============================================")
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestApp)
    RET = unittest.TextTestRunner(verbosity=2).run(SUITE)
    # print (out the test results)
    print (str(RET))
    if RET.printErrors() != None:
        print ("Error: " + RET.printErrors())

    if RET.wasSuccessful():
        exit(0)
    else:
        exit(1)
