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

        print "Results are:", results

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
            print "   Result for", i, "is", result_counts[i]

        for i in range(len(result_counts)):
            # Ensure results are within 10% of expected
            self.assertLess(abs(expected_result_counts[i] - result_counts[i]),
                            loop_count/10,
                            msg="Results out of range for instance " + str(i) +
                            " Actual result: " + str(result_counts[i]) +
                            " Expected result: " + str(expected_result_counts[i]))

    # End of test_RandomWeighted()

    def test_GenProfileFromMotifs(self):

        motifs = [
                ['T', 'A', 'A', 'C'],
                ['G', 'T', 'C', 'T'],
                ['A', 'C', 'T', 'A'],
                ['A', 'G', 'G', 'T']
            ]

        profile = utils.GenProfileFromMotifs(motifs)

        print "Profile returned is:"
        for row in profile:
            print "   ", row

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

        print "Found:", probabilities

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

        print "  Kmer:   Results:  Expected Results:"
        for kmer in utils.KmersFromText(text, k):
            print "   {kmer}   {result:6}       {expected_result:6}".format(
                  kmer=kmer, result=results[kmer],
                  expected_result=int(expected_results[kmer]))

        for kmer in sorted(expected_results.keys()):
            self.assertLess(abs(results[kmer] - int(expected_results[kmer])),
                            iter_count/50,
                            msg="Results out of range for '" + kmer +
                            "'. Actual result: " + str(results[kmer]) +
                            ", Expected result: " + str(expected_results[kmer]))

    # End of test_ProfileRandomlyGeneratedKmer()

# End of class TestApp


if __name__ == "__main__":

    LOG.debug("Unittest args: " + ' '.join(sys.argv))
    LOG.debug("Begining of the test suite:")
    LOG.debug("===============================================")
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestApp)
    RET = unittest.TextTestRunner(verbosity=2).run(SUITE)
    # Print out the test results
    print str(RET)
    if RET.printErrors() != None:
        print "Error: " + RET.printErrors()

    if RET.wasSuccessful():
        exit(0)
    else:
        exit(1)
