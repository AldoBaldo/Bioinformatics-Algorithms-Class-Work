#!/usr/bin/env python

import sys
import argparse
import array
import itertools
import time

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

def PatternCount(sequence, pattern, max_hamming_distance=0):

    count = 0
    len_sequence = len(sequence)
    len_pattern = len(pattern)
    srch_count = len_sequence - len_pattern + 1

    for i in range(0, srch_count):
        if FindHammingDistance(pattern, sequence[i:i+len_pattern]) <= max_hamming_distance:
            count += 1

    return count

# End of PatternCount()

# This array converts the lowest three bits of a base letter (A,C,G,T) to number that can be used
# to index into an array.
#     'A' is 0x41, so the lowest 3 bits are 001, so index 1 maps to a value of 0.
#     'C' is 0x43, so the lowest 3 bits are 011, so index 3 maps to a value of 1.
#     'G' is 0x47, so the lowest 3 bits are 111, so index 7 maps to a value of 2.
#     'T' is 0x54, so the lowest 3 bits are 100, so index 4 maps to a value of 3.

BaseToNumberArray = [0, 0, 0, 1, 3, 0, 0, 2]

def BaseToNumber(base):
    if base not in ['A', 'C', 'G', 'T']:
        raise Exception("Invalid base")
    return BaseToNumberArray[ord(base) & 0b111]
# End of BaseToNumber()

def PatternToNumber(pattern):

    result = 0

    for c in pattern:
        result = (result<<2) + BaseToNumber(c)

    return result

# End of PatternToNumber()

NumberToBaseArray           = ['A', 'C', 'G', 'T']
NumberToBaseComplimentArray = ['T', 'G', 'C', 'A']

def NumberToPattern(number, k):

    result = ''
    for i in range(2*k-2, -2, -2):
        window = 0b11 << i
        result += NumberToBaseArray[(number & window)>>i]

    return result

# End of NumberToPattern()

def FindMostFrequentUsingDict(sequence, k, t=0):

    '''
    If t is non-zero, find all patterns with frequency greater than t.
    Otherwise, find the pattern that occurs the most.
    '''

    pattern_counts = {}
    len_sequence = len(sequence)
    srch_count = len_sequence - k + 1

    for i in range(0, srch_count):
        cur_pattern = sequence[i:i+k]
        if cur_pattern in pattern_counts:
            pattern_counts[cur_pattern] += 1
        else:
            pattern_counts[cur_pattern] = 1

    max_count = t
    result = []
    for pattern, count in pattern_counts.iteritems():
        if count > max_count:
            if t == 0:
                max_count = count
            result = [pattern]
        elif count == max_count:
            if pattern not in result:
                result.append(pattern)

    return sorted(result)

# End of FindMostFrequentUsingDict()

def FindMostFrequentUsingArray(sequence, k, t=0):

    '''
    If t is non-zero, find all patterns with frequency greater than t.
    Otherwise, find the pattern that occurs the most.
    '''

    pattern_counts = array.array('I', itertools.repeat(0, 4**k))
    len_sequence = len(sequence)
    srch_count = len_sequence - k + 1

    for i in range(0, srch_count):
        cur_pattern = sequence[i:i+k]
        pattern_counts[PatternToNumber(cur_pattern)] += 1

    max_count = t
    result = []
    for index, count in enumerate(pattern_counts):
        if count > max_count:
            if t > 0:
                max_count = count
            result = []    # Remove old results
        if count == max_count:
            pattern = NumberToPattern(index, k)
            if pattern not in result:
                result.append(pattern)

    return sorted(result)

# End of FindMostFrequentUsingArray()

def ReverseCompliment(pattern):

    result = ''
    for ch in pattern:
        result += NumberToBaseComplimentArray[BaseToNumber(ch)]
    result = result[::-1]   # Reverse the string
    return result

# End of ReverseCompliment()

def FindPattern(pattern, genome):

    cur_index = 0
    result = []
    while True:
        cur_index = genome.find(pattern, cur_index)
        if cur_index >= 0:
            result.append(cur_index)
            cur_index += 1
        else:
            break

    return result

# End of FindPattern()

def FindPatternClump(k, L, t, genome):

    window_start = 0
    window_end = L-1
    search_end = len(genome) - k
    patterns_found = []

    while window_end < search_end:
        new_patterns = FindMostFrequentUsingDict(genome[window_start:window_end], k, t)
        for pattern in new_patterns:
            if pattern not in patterns_found:
                patterns_found.append(pattern)

        window_start += 1
        window_end += 1

    return patterns_found

# End of FindPatternClump()

def FindMinimumSkews(genome):

    cur_skew = 0
    min_skew = 0
    min_skew_locci = []

    debug_skews = []

    for i, bp in enumerate(genome):
        if bp == 'C':
            cur_skew -= 1
        elif bp == 'G':
            cur_skew += 1

        if cur_skew == min_skew:
            min_skew_locci.append(i+1)
        elif cur_skew < min_skew:
            min_skew_locci = [i+1]
            min_skew = cur_skew

        debug_skews.append(cur_skew)

    # print '  ' + ',   '.join(genome)
    # print ', '.join(['% 3d' % x for x in debug_skews])
    # print ', '.join(['% 3d' % x for x in range(1, len(genome)+1)])

    return min_skew_locci

# End of FindMinimumSkews()

def FindHammingDistance(pattern_a, pattern_b):

    hamming_distance = 0

    for i in range(0, len(pattern_a)):
        if pattern_a[i] != pattern_b[i]:
            hamming_distance += 1

    return hamming_distance

# End of FindHammingDistance()

def FindSumOfHammingDistances(pattern, dna):

    hamming_sum = 0
    k = len(pattern)
    n = len(dna[0])    # Assume all strands in dna are the same length

    for strand in dna:
        min_hamming = k
        for i in range(n-k+1):
            motif = strand[i:i+k]
            hamming_distance = FindHammingDistance(pattern, motif)
            min_hamming = min(min_hamming, hamming_distance)
            if min_hamming == 0:
                break
        hamming_sum += min_hamming
        print "  For pattern <" + pattern + "> in strand <" + strand + ">, the min hamming distance is:", min_hamming

    return hamming_sum

# End of FindSumOfHammingDistances()

def FindApproximatePatternMatches(pattern, max_hamming_distance, text):

    result = []
    pattern_len = len(pattern)

    for i in range(0, len(text) - pattern_len + 1):

        cur_hamming_distance = FindHammingDistance(pattern, text[i:i+pattern_len])
        if cur_hamming_distance <= max_hamming_distance:
            result.append(i)

    return result

# End of FindApproximatePatternMatches()

def FindMostFrequentWithMismatches(k, d, text):

    '''Look for all possible patterns in the given text'''

    pattern_counts = []
    max_count_found = 0

    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        pattern_count = PatternCount(text, pattern, d)
        if pattern_count == max_count_found:
            pattern_counts.append(pattern)
        elif pattern_count > max_count_found:
            pattern_counts = [pattern]
            max_count_found = pattern_count

    return pattern_counts

# End of FindMostFrequentWithMismatches()

def FindMostFrequentWithMismatchesAndReverseCompiliment(k, d, text):

    '''Look for all possible patterns in the given text'''

    pattern_counts = []
    max_count_found = 0

    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        reverse_compliment = ReverseCompliment(pattern)
        pattern_count = PatternCount(text, pattern, d)
        pattern_count += PatternCount(text, reverse_compliment, d)
        if pattern_count == max_count_found:
            pattern_counts.append(pattern)
        elif pattern_count > max_count_found:
            pattern_counts = [pattern]
            max_count_found = pattern_count

    return pattern_counts

# End of FindMostFrequentWithMismatchesAndReverseCompiliment()


FindOtherBasesArray = [
    None,
    ['C', 'G', 'T'],   # All bases other than 'A' (lowest 3 bits are 001)
    None,
    ['A', 'G', 'T'],   # All bases other than 'C' (lowest 3 bits are 011)
    ['A', 'C', 'G'],   # All bases other than 'T' (lowest 3 bits are 100)
    None,
    None,
    ['A', 'C', 'T']    # All bases other than 'G' (lowest 3 bits are 111)
]
def FindOtherBases(base):

    return FindOtherBasesArray[ord(base) & 0b111]

# End of FindOtherBases

def FindImmediateNeighbors(pattern):

    neighborhood = [pattern]

    for i in range(len(pattern)):
        for base in FindOtherBases(pattern[i]):
            neighborhood.append(pattern[:i] + base + pattern[i+1:])

    return neighborhood

# End of FindImmediateNeighbors()

def FindNeighbors(pattern, d, next_base_list=['A', 'C', 'G', 'T']):

    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []
    suffix_pattern = pattern[1:]
    suffix_neighbors = FindNeighbors(suffix_pattern, d, next_base_list=next_base_list)
    for suf_neighbor in suffix_neighbors:
        if FindHammingDistance(suffix_pattern, suf_neighbor) < d:
            for base in next_base_list:
                neighborhood.append(base + suf_neighbor)
        else:
            neighborhood.append(pattern[0] + suf_neighbor)

    return neighborhood

# End of FindNeighbors()

def FindPatternFrequencies(k, text):

    pattern_counts = array.array('I', itertools.repeat(0, 4**k))
    len_text = len(text)
    srch_count = len_text - k + 1

    for i in range(0, srch_count):
        cur_pattern = text[i:i+k]
        pattern_counts[PatternToNumber(cur_pattern)] += 1

    return pattern_counts

# End of FindPatternFrequencies()

def FindImplantedMotifs(k, d, dna):

    pattern_lists = []

    # Find all of the patterns per dna strand
    for strand in dna:
        pattern_list = []
        for i in range(len(strand)-k+1):
            pattern = strand[i:i+k]
            pattern_list += FindNeighbors(pattern, d)
        pattern_lists.append(pattern_list)

    # Find patterns that are in all lists
    first_list = pattern_lists[0]
    remaining_lists = pattern_lists[1:]
    unique_patterns = []
    for pattern in first_list:
        pattern_in_all_lists = True
        for other_list in remaining_lists:
            if pattern not in other_list:
                pattern_in_all_lists = False
        if pattern_in_all_lists and pattern not in unique_patterns:
            unique_patterns.append(pattern)

    return sorted(unique_patterns)

# End of FindImplantedMotifs()

def FindMedianString(k, dna):

    t = len(dna)
    distance = k * t
    median = ''
    saved_median_distances = {}
    saved_medians_of_min_distance = []

    # Find distance for each possible k-mer
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        print "Finding hamming sum for <" + pattern + ">"
        new_distance = FindSumOfHammingDistances(pattern, dna)
        print "New hamming sum is:", new_distance
        saved_median_distances[pattern] = new_distance
        if new_distance < distance:
            distance = new_distance
            median = pattern
            saved_medians_of_min_distance = [pattern]
        elif new_distance == distance:
            saved_medians_of_min_distance.append(pattern)

    for pattern in sorted(saved_median_distances):
        print "The hamming sum for >" + pattern + "> is:", saved_median_distances[pattern]

    print "The following patterns share a min hamming distance of:", distance, ":", ' '.join(saved_medians_of_min_distance)

    return median

# End of FindMedianString()

def ProfileMostProbableKmer(text, k, profile):

    result = ''
    best_probability = 0

    for i in range (len(text) - k + 1):
        kmer = text[i:i+k]
        cur_probability = reduce(lambda x, y: x*y,
                                 [profile[BaseToNumber(base)][i] for i, base in enumerate(kmer)])

        if cur_probability > best_probability:
            best_probability = cur_probability
            result = kmer

    return result

# End of ProfileMostProbableKmer()

def GenProfileFromMotifs(motifs):

    num_motifs = len(motifs)
    motif_len = len(motifs[0])
    profile = list(itertools.repeat(list(itertools.repeat(0.0, motif_len)), 4))

    for i in range(num_motifs):
        for j in range(motif_len):
            profile[i, BaseToNumber(motifs[i,j])] += 1.0/num_motifs

    return profile

# End of GenProfileFromMotifs()

def CalculateScore(motifs):

    num_motifs = len(motifs)
    motif_len = len(motifs[0])
    score = 0

    for i in range(motif_len):
        counts = [0,0,0,0]
        for j in range(num_motifs):
            base = motifs[j][i]
            counts[BaseToNumber(base)] += 1
        score += max(counts)

    return score

# End of CalculateScore()

def GreedyMotifSearch(k, t, dna):

    best_motifs = [x[:k] for x in dna]

    for k in range(len(dna[0]) - k + 1):
        motifs = [dna[0][k:k+1]]

        for i in range(1:t):

            profile = GenProfileFromMotifs(motifs)
            motifs.append(ProfileMostProbableKmer(dna[i], k, profile))

        if CalculateScore(motifs) < CalculateScore(best_motifs):
            best_motifs = motifs

    return best_motifs

# End of GreedyMotifSearch()

def Exercise_ba1b_FindMostFrequentString():

    parser = argparse.ArgumentParser(description="Find the most frequent patterns of length k")
    parser.add_argument('k', type=int, help="The length of the pattern to search for")
    parser.add_argument('sequence', type=str, help="The DNA sequence to search")

    args = parser.parse_args()

    start_time = time.time()
    result = FindMostFrequentUsingDict(args.sequence, args.k)
    end_time = time.time()
    print "Dict Result is:\n", ' '.join(result)
    print "Elapsed time = ", end_time - start_time

    start_time = time.time()
    result = FindMostFrequentUsingArray(args.sequence, args.k)
    end_time = time.time()
    print "Array Result is:\n", ' '.join(result)
    print "Elapsed time = ", end_time - start_time

# End of Exercise_ba1b_FindMostFrequentString()

def Exercise_ba1c_ReverseCompliment():

    parser = argparse.ArgumentParser(description="Find the reverse compliment of a pattern")
    parser.add_argument('sequence', type=str, help="The DNA sequence to search")

    args = parser.parse_args()

    result = ReverseCompliment(args.sequence)
    print 'The reverse compliment of "' + args.sequence + '" is "' + result + '"'

# End of Exercise_ba1c_ReverseCompliment()

def Exercise_ba1d_FindPatternInGenome():

    parser = argparse.ArgumentParser(description="Find the most frequent patterns of length k")
    parser.add_argument('pattern', type=str, help="The DNA pattern to search for")
    parser.add_argument('genome', type=str, help="The genome to search for the pattern in")

    args = parser.parse_args()

    result = FindPattern(args.pattern, args.genome)
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

    result = FindPatternClump(args.k, args.L, args.t, args.genome)
    print 'The distinct k-mers forming (L,t)-clumps are:', ' '.join([str(x) for x in result])

# End of Exercise_ba1d_FindPatternClumpInGenome()

def Exercise_ba1f_FindMinimumSkews():
    '''Find Find the minimum GC skews in the given genome.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1e_FindMinimumSkews.__doc__)
    parser.add_argument('genome', type=str, help="The genome to find minimum skews in")

    args = parser.parse_args()

    result = FindMinimumSkews(args.genome)
    print 'The minimum skews for this genome are:', ' '.join([str(x) for x in result])

# End of Exercise_ba1e_FindMinimumSkews()

def Exercise_ba1g_FindHammingDistance():
    '''Find the hamming distance between two sequences.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1g_FindHammingDistance.__doc__)
    parser.add_argument('pattern_a', type=str, help="The first sequence")
    parser.add_argument('pattern_b', type=str, help="The second sequence")

    args = parser.parse_args()

    result = FindHammingDistance(args.pattern_a, args.pattern_b)
    print 'The hamming distance between these two patters is:', result

# End of Exercise_ba1g_FindHammingDistance()

def Exercise_ba1h_FindApproximatePatternMatches():
    '''Find all instances of a pattern in a text with a Hamming distance <= d.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1h_FindApproximatePatternMatches.__doc__)
    parser.add_argument('pattern', type=str, help="The pattern to look for")
    parser.add_argument('max_hamming_distance', type=int, help="The max Hamming distance")
    parser.add_argument('text', type=str, help="The text to search in")

    args = parser.parse_args()

    result = FindApproximatePatternMatches(args.pattern, args.max_hamming_distance, args.text)
    print 'Approximate matches of the pattern can be found at:', ' '.join([str(x) for x in result])

# End of Exercise_ba1h_FindApproximatePatternMatches()

def Exercise_ba1i_FindMostFrequentWithMismatches():
    '''Find the most frequent k-mers in a text with at most d mismatches.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1i_FindMostFrequentWithMismatches.__doc__)
    parser.add_argument('k', type=int, help="The length of the sequence to search for")
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('text', type=str, help="The text to search in")

    args = parser.parse_args()

    result = FindMostFrequentWithMismatches(args.k, args.d, args.text)
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

    result = FindMostFrequentWithMismatchesAndReverseCompiliment(args.k, args.d, args.text)
    print 'The most frequent k-mers in the text with mismatches and reverse compliments are:', ' '.join(result)

# End of Exercise_ba1j_FindMostFrequentWithMismatchesAndReverseCompiliment()

def Exercise_ba1k_CountingFrequencies():
    '''Return a frequency array with the frequencies of all patterns of length k.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1k_CountingFrequencies.__doc__)
    parser.add_argument('k', type=int, help="The length of the pattern to count")
    parser.add_argument('text', type=str, help="The text to find the patterns in")

    args = parser.parse_args()

    result = FindPatternFrequencies(args.k, args.text)
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

    result = PatternToNumber(args.pattern)
    print "The unique number for the pattern is:", result

    # Test the result
    pattern_back = NumberToPattern(result, len(args.pattern))
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

    result = NumberToPattern(args.n, args.k)
    print 'The pattern for n, given k, is:', result

    # Test the result
    number_back = PatternToNumber(result)
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

    result = FindNeighbors(args.pattern, args.d)
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

    result = FindImplantedMotifs(args.k, args.d, args.dna)
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

    result = FindMedianString(args.k, args.dna)
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

    result = ProfileMostProbableKmer(text, k, profile)
    print 'The profile-most probable k-mer is:', result

# End of Exercise_ba2c_ProfileMostProbableKmer()

def Exercise_ba2d_GreedyMotifSearch():
    '''Find the most probably k-mer in text for the given profile.'''

    print "Enter the data (text, k, profile matrix):"

    k = int(raw_input())
    t = int(raw_input())
    dna = []
    while True:
        line = raw_input()
        if len(line) == 0:
            break
        else:
            dna.append([float(x) for x in line.split()])

    result = GreedyMotifSearch(k, t, dna)
    print 'The best motifs for the given dna is:'
    for row in result:
        print '   ' + row

# End of Exercise_ba2d_GreedyMotifSearch()

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
        result = FindNeighbors("ACG", 1, next_base_list = combo)
        print "For combo:", combo
        print "                Expected results are:", expected_results
        print "                My results are      :", result
        if result == expected_results:
            winners.append(combo)

    print "The winning combos are:", winners

# End of TestCombo()


if __name__ == "__main__":

    # Test Code

    #import pdb; pdb.set_trace()
    #PatternToNumber('GT')
    #NumberToPattern(11,2)
    #NumberToPattern(PatternToNumber('GT'), 2)

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
    Exercise_ba2d_GreedyMotifSearch()
