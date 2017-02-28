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
