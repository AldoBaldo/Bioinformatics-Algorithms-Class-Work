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
