#!/usr/bin/env python

import sys
import argparse
import array
import itertools
import time

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

def FindNeighbors(pattern, d):

    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []
    suffix_pattern = pattern[1:]
    suffix_neighbors = FindNeighbors(suffix_pattern, d)
    for suf_neighbor in suffix_neighbors:
        if FindHammingDistance(suffix_pattern, suf_neighbor) < d:
            for base in ['A', 'C', 'G', 'T']:
            #for base in ['T', 'G', 'C', 'A']:
            #for base in ['T', 'G', 'C', 'A']:
                neighborhood.append(base + suf_neighbor)
        else:
            neighborhood.append(pattern[0] + suf_neighbor)

    return neighborhood

# End of FindNeighbors()

def FindMostFrequentWithMismatchesBetter(k, d, text):

    pass

# End of FindMostFrequentWithMismatchesFaster()

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

    result = FindMostFrequentWithMismatchesBetter(args.k, args.d, args.text)
    print 'The most frequent k-mers in the text with mismatches are:', ' '.join(result)

# End of Exercise_ba1i_FindMostFrequentWithMismatches()

def Exercise_ba1n_FindNeighbors():
    '''Find the most frequent k-mers in a text with at most d mismatches.'''

    parser = argparse.ArgumentParser(description=Exercise_ba1n_FindNeighbors.__doc__)
    parser.add_argument('d', type=int, help="The max Hamming distance")
    parser.add_argument('pattern', type=str, help="The pattern to find neighbors for")

    args = parser.parse_args()

    result = FindNeighbors(args.pattern, args.d)
    print 'The neighbors of', args.pattern, 'are:', ' '.join(result)

# End of Exercise_ba1n_FindNeighbors()



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

    Exercise_ba1n_FindNeighbors()