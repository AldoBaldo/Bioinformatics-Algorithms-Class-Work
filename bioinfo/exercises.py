#!/usr/bin/env python

import sys
import argparse
import time

import utils

def TestResults(my_results, expected_results):

    error_count = 0

    if len(my_results) != len(expected_results):
        print "Results differ in length.", len(my_results), 'vs', len(expected_results)
        return

    for item in my_results:
        if item not in expected_results:
            error_count += 1
            print "Item", item, "from my results is not found in expected results"

    for item in expected_results:
        if item not in my_results:
            error_count += 1
            print "Item", item, "from expected results is not found in my results"

    print error_count, "errors found"

# End of Test Results()


def FindUniqueBaseLists(start = ['A', 'C', 'G', 'T']):

    if len(start) <= 1:
        return [start]

    result = []

    for base in start:
        new_pool = start[:]
        new_pool.remove(base)
        result += [[base] + x for x in FindUniqueBaseLists(new_pool)]

    return result

# End of FindUniqueBaseLists()

def TestCombo():

    combos = FindUniqueBaseLists()
    expected_results = ["CCG", "TCG", "GCG", "AAG", "ATG", "AGG", "ACA", "ACC", "ACT", "ACG"]
    winners = []

    for combo in combos:
        result = utils.FindNeighbors("ACG", 1, next_base_list = combo)
        print "For combo:", combo
        print "                Expected results are:", expected_results
        print "                My results are      :", result
        if result == expected_results:
            winners.append(combo)

    print "The winning combos are:", winners

# End of TestCombo()

def Exercise_ba1a_PatternCount():

    '''Count how often the given pattern appears in the given text.'''

    print "Enter the data (text, pattern):"

    text = raw_input()
    pattern = raw_input()

    print "Python version:", sys.version
    print "Pattern Length:", len(pattern)
    print "Text Length:", len(text)

    result = utils.PatternCount(text, pattern)

    print 'The pattern is repeated the following number of times:'
    print '   ', result

    print "Enter the results to compare against:"
    expected_result = raw_input()
    if len(expected_result) > 0:
        expected_result = int(expected_result)
        if expected_result != result:
            print "Results don't match.  Expected", expected_result, "but found", result
        else:
            print "Results match"
    else:
        print "Skipping results validation."

# End of Exercise_ba1a_PatternCount()

def Exercise_ba1b_FindMostFrequentString():

    parser = argparse.ArgumentParser(description="Find the most frequent patterns of length k")
    parser.add_argument('k', type=int, help="The length of the pattern to search for")
    parser.add_argument('sequence', type=str, help="The DNA sequence to search")

    args = parser.parse_args()

    start_time = time.time()
    result = utils.FindMostFrequentUsingDict(args.sequence, args.k)
    end_time = time.time()
    print "Dict Result is:\n", ' '.join(result)
    print "Elapsed time = ", end_time - start_time

    start_time = time.time()
    result = utils.FindMostFrequentUsingArray(args.sequence, args.k)
    end_time = time.time()
    print "Array Result is:\n", ' '.join(result)
    print "Elapsed time = ", end_time - start_time

# End of Exercise_ba1b_FindMostFrequentString()

def Exercise_ba1c_ReverseCompliment():

    parser = argparse.ArgumentParser(description="Find the reverse compliment of a pattern")
    parser.add_argument('sequence', type=str, help="The DNA sequence to search")

    args = parser.parse_args()

    result = utils.ReverseCompliment(args.sequence)
    print 'The reverse compliment of "' + args.sequence + '" is "' + result + '"'

# End of Exercise_ba1c_ReverseCompliment()

def Exercise_ba1d_FindPatternInGenome():

    parser = argparse.ArgumentParser(description="Find the most frequent patterns of length k")
    parser.add_argument('pattern', type=str, help="The DNA pattern to search for")
    parser.add_argument('genome', type=str, help="The genome to search for the pattern in")

    args = parser.parse_args()

    result = utils.FindPattern(args.pattern, args.genome)
    print 'The pattern was found at the following locations:', ' '.join([str(x) for x in result])

# End of Exercise_ba1d_FindPatternInGenome()

def Exercise_ba1e_FindPatternClumpInGenome():
    '''Find (L, t)-clumps of k-mers in a genome'''

    parser = argparse.ArgumentParser(description=Exercise_ba1d_FindPatternClumpInGenome.__doc__)
    parser.add_argument('k', type=int, help="The length of the pattern to look for")
    parser.add_argument('L', type=int, help="The size of the sliding window")
    parser.add_argument('t', type=int, help="The minimum number of copies for a pattern to repeat")
    parser.add_argument('genome', type=str, help="The genome to search for the pattern in")

    args = parser.parse_args()

    result = utils.FindPatternClump(args.k, args.L, args.t, args.genome)
    print 'The distinct k-mers forming (L,t)-clumps are:', ' '.join([str(x) for x in result])

# End of Exercise_ba1d_FindPatternClumpInGenome()

def Exercise_ba1f_FindMinimumSkews():
    '''Find Find the minimum GC skews in the given genome.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1e_FindMinimumSkews.__doc__)
    parser.add_argument('genome', type=str, help="The genome to find minimum skews in")

    args = parser.parse_args()

    result = utils.FindMinimumSkews(args.genome)
    print 'The minimum skews for this genome are:', ' '.join([str(x) for x in result])

# End of Exercise_ba1e_FindMinimumSkews()

def Exercise_ba1g_FindHammingDistance():
    '''Find the hamming distance between two sequences.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1g_FindHammingDistance.__doc__)
    parser.add_argument('pattern_a', type=str, help="The first sequence")
    parser.add_argument('pattern_b', type=str, help="The second sequence")

    args = parser.parse_args()

    result = utils.FindHammingDistance(args.pattern_a, args.pattern_b)
    print 'The hamming distance between these two patters is:', result

# End of Exercise_ba1g_FindHammingDistance()

def Exercise_ba1h_FindApproximatePatternMatches():
    '''Find all instances of a pattern in a text with a Hamming distance <= d.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1h_FindApproximatePatternMatches.__doc__)
    parser.add_argument('pattern', type=str, help="The pattern to look for")
    parser.add_argument('max_hamming_distance', type=int, help="The max Hamming distance")
    parser.add_argument('text', type=str, help="The text to search in")

    args = parser.parse_args()

    result = utils.FindApproximatePatternMatches(args.pattern, args.max_hamming_distance, args.text)
    print 'Approximate matches of the pattern can be found at:', ' '.join([str(x) for x in result])

# End of Exercise_ba1h_FindApproximatePatternMatches()

def Exercise_ba1i_FindMostFrequentWithMismatches():
    '''Find the most frequent k-mers in a text with at most d mismatches.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1i_FindMostFrequentWithMismatches.__doc__)
    parser.add_argument('k', type=int, help="The length of the sequence to search for")
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('text', type=str, help="The text to search in")

    args = parser.parse_args()

    result = utils.FindMostFrequentWithMismatches(args.k, args.d, args.text)
    print 'The most frequent k-mers in the text with mismatches are:', ' '.join(result)

# End of Exercise_ba1i_FindMostFrequentWithMismatches()

def Exercise_ba1j_FindMostFrequentWithMismatchesAndReverseCompiliment():
    '''Find the most frequent k-mers in a text with at most d mismatches,
    and including their reverse compliments.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1j_FindMostFrequentWithMismatchesAndReverseCompiliment.__doc__)
    parser.add_argument('k', type=int, help="The length of the sequence to search for")
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('text', type=str, help="The text to search in")

    args = parser.parse_args()

    result = utils.FindMostFrequentWithMismatchesAndReverseCompiliment(args.k, args.d, args.text)
    print 'The most frequent k-mers in the text with mismatches and reverse compliments are:', ' '.join(result)

# End of Exercise_ba1j_FindMostFrequentWithMismatchesAndReverseCompiliment()

def Exercise_ba1k_CountingFrequencies():
    '''Return a frequency array with the frequencies of all patterns of length k.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1k_CountingFrequencies.__doc__)
    parser.add_argument('k', type=int, help="The length of the pattern to count")
    parser.add_argument('text', type=str, help="The text to find the patterns in")

    args = parser.parse_args()

    result = utils.FindPatternFrequencies(args.k, args.text)
    print 'The frequency array for the given k and text are:', ' '.join([str(x) for x in result])

    # Test the results when the results are known
#    with open('ba1k_results.txt', 'r') as expected_results_fh:
#        expected_results = []
#        for line in expected_results_fh.readlines():
#            expected_results += [int(x) for x in line.split()]
#
#    TestResults(result, expected_results)

# End of Exercise_ba1k_CountingFrequencies()

def Exercise_ba1l_PatternToNumber():

    '''Return a unique number for the specified pattern.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1k_CountingFrequencies.__doc__)
    parser.add_argument('pattern', type=str, help="The pattern to index")

    args = parser.parse_args()

    result = utils.PatternToNumber(args.pattern)
    print "The unique number for the pattern is:", result

    # Test the result
    pattern_back = utils.NumberToPattern(result, len(args.pattern))
    if pattern_back == args.pattern:
        print "Test passed for PatternToNumber"
    else:
        print "Test failed for PatternToNumber"

# End of Exercise_ba1l_PatternToNumber()

def Exercise_ba1m_NumberToPattern():

    '''Reverse the "PatternToNumber" transformation and retrieve the correct pattern
    for the given number, given k.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1k_CountingFrequencies.__doc__)
    parser.add_argument('n', type=int, help="The number for the pattern")
    parser.add_argument('k', type=int, help="The length of the sequence used in calculating n.")

    args = parser.parse_args()

    result = utils.NumberToPattern(args.n, args.k)
    print 'The pattern for n, given k, is:', result

    # Test the result
    number_back = utils.PatternToNumber(result)
    if number_back == args.n and len(result) == args.k:
        print "Test passed for NumberToPattern"
    else:
        print "Test failed for NumberToPattern"

# End of Exercise_ba1m_NumberToPattern()

def Exercise_ba1n_FindNeighbors():
    '''Find the most frequent k-mers in a text with at most d mismatches.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1n_FindNeighbors.__doc__)
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('pattern', type=str, help="The pattern to find neighbors for")

    args = parser.parse_args()

    result = utils.FindNeighbors(args.pattern, args.d)
    print 'The neighbors of', args.pattern, 'are:', ' '.join(result)

    #expected_results = ["CCG", "TCG", "GCG", "AAG", "ATG", "AGG", "ACA", "ACC", "ACT", "ACG"]

    #with open('ba1n_test_output.txt', 'r') as expected_results_fh:
    #    expected_results = [x.strip() for x in expected_results_fh.readlines() if len(x.strip()) > 0]

    #TestResults(result, expected_results)

# End of Exercise_ba1n_FindNeighbors()

def Exercise_ba2a_FindImplantedMotifs():
    '''Find k-mer motifs with at most d mismatches in the given DNA sequences.'''

    parser = argparse.ArgumentParser(description=Exercise_ba2a_FindImplantedMotifs.__doc__)
    parser.add_argument('k', type=int, help="The length of patterns to look for")
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('dna', type=str, nargs='+', help="A list of DNA sequences to search in")

    args = parser.parse_args()

    result = utils.FindImplantedMotifs(args.k, args.d, args.dna)
    print 'The implanted motifs are', ' '.join(result)

#    with open('ba2a_2_test_results.txt', 'r') as expected_results_fh:
#        expected_results = expected_results_fh.readline().split()
#    TestResults(result, expected_results)

# End of Exercise_ba2a_FindImplantedMotifs()

def Exercise_ba2b_FindMedianString():
    '''Find a string that minimizes the hamming distance between itself and any
    pattern in a list of DNA sequences.'''

    parser = argparse.ArgumentParser(description=Exercise_ba2b_FindMedianString.__doc__)
    parser.add_argument('k', type=int, help="The length of patterns to look for")
    parser.add_argument('dna', type=str, nargs='+', help="A list of DNA sequences to search in")

    args = parser.parse_args()

    result = utils.FindMedianString(args.k, args.dna)
    print 'The median string is:', result

# End of Exercise_ba2b_FindMedianString()

def Exercise_ba2c_ProfileMostProbableKmer():
    '''Find the most probably k-mer in text for the given profile.'''

    print "Enter the data (text, k, profile matrix):"

    text = raw_input()
    k = int(raw_input())
    profile = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            profile.append([float(x) for x in line.split()])

    result = utils.ProfileMostProbableKmer(text, k, profile)
    print 'The profile-most probable k-mer is:', result

# End of Exercise_ba2c_ProfileMostProbableKmer()

def Exercise_ba2e_GreedyMotifSearch():
    '''Find the best k-mers motifs from the given DNA using the
    GreedyMotifSearch method.'''

    print "Enter the data (k t, dna):"

    k, t = [int(x) for x in raw_input().split()]
    # t = int(raw_input())
    dna = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            dna.append(line)

    result = utils.GreedyMotifSearch(k, t, dna)
    print 'The best motifs for the given dna is:'
    for row in result:
        print ('   ' + row)

    print "Enter the results to compare against:"
    expected_results = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            expected_results.append(line)
    if len(expected_results) > 0:
        TestResults(result, expected_results)
    else:
        print "Skipping results validation."

# End of Exercise_ba2e_GreedyMotifSearch()

def Exercise_ba2f_RandomizedMotifSearch():
    '''Find the best k-mers motifs from the given DNA using the
    RandomizedMotifSearch method.'''

    print "Enter the data (k t, dna):"

    k, t = [int(x) for x in raw_input().split()]
    # t = int(raw_input())
    dna = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            dna.append(line)

    best_result = None
    for i in range(1000):
        result = utils.RandomizedMotifSearch(k, t, dna)
        if (best_result == None
                or CalculateScore(result) < CalculateScore(best_result)):
            best_result = result
        #debug:
        if i%10 == 0:
            # print 'Round', i, best_result
            print '.',

    print 'The best motifs for the given dna is:'
    for row in best_result:
        print '   ' + row

    print "Enter the results to compare against:"
    expected_results = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            expected_results.append(line)
    if len(expected_results) > 0:
        TestResults(best_result, expected_results)
    else:
        print "Skipping results validation."

# End of Exercise_ba2f_RandomizedMotifSearch()

def Exercise_ba2g_GibbsSampler():
    '''Find the best k-mers motifs from the given DNA using the
    Gibbs Sampler method.'''

    print "Enter the data (k t N, dna):"

    k, t, N = [int(x) for x in raw_input().split()]
    # t = int(raw_input())
    dna = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            dna.append(line)

    result = utils.GibbsSampler(dna, k, t, N)

    print 'The best motifs for the given dna is:'
    for row in result:
        print '   ' + row

    print ("Enter the results to compare against:")
    expected_results = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            expected_results.append(line)
    if len(expected_results) > 0:
        TestResults(result, expected_results)
    else:
        print "Skipping results validation."

# End of Exercise_ba2g_GibbsSampler()


if __name__ == "__main__":

    # Test Code

    #import pdb; pdb.set_trace()
    #PatternToNumber('GT')
    #NumberToPattern(11,2)
    #NumberToPattern(PatternToNumber('GT'), 2)

    Exercise_ba1a_PatternCount()
    # Exercise_ba1b_FindMostFrequentString()
    # Exercise_ba1c_ReverseCompliment()
    # Exercise_ba1d_FindPatternInGenome()
    # Exercise_ba1e_FindPatternClumpInGenome()
    # Exercise_ba1f_FindMinimumSkews()
    # Exercise_ba1g_FindHammingDistance()
    # Exercise_ba1h_FindApproximatePatternMatches()
    # Exercise_ba1i_FindMostFrequentWithMismatches()
    # Exercise_ba1j_FindMostFrequentWithMismatchesAndReverseCompiliment()
    # Exercise_ba1k_CountingFrequencies()
    # Exercise_ba1l_PatternToNumber()
    # Exercise_ba1m_NumberToPattern()
    # Exercise_ba1n_FindNeighbors()

    # Exercise_ba2a_FindImplantedMotifs()
    # Exercise_ba2b_FindMedianString()
    # Exercise_ba2c_ProfileMostProbableKmer()
    # Exercise_ba2e_GreedyMotifSearch()
    # Exercise_ba2f_RandomizedMotifSearch()
    # Exercise_ba2g_GibbsSampler()
