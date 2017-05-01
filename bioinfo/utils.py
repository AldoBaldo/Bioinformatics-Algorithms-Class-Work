#!/usr/bin/env python

import sys
import argparse
import array
import itertools
import random
import time
from collections import OrderedDict
from functools import reduce
import copy

DEBUG = False

def debug(*args):
    if DEBUG:
        print(*args)

def GetInput(*args):

    debug("GetInput(" + str(args) + ")",)

    result = []
    return_list = (len(args) > 1)

    for arg in args:
        if arg in ['k', 't', 'd', 'i']:
            result.append(int(input().strip()))
        elif arg in ['k t', 'k d', 'i i']:
            return_list = True
            result += [int(x.strip()) for x in input().split()]
        elif arg in ['pattern', 's']:
            result.append(input().strip())
        elif arg in ['s s']:
            return_list = True
            result.append(input().split())
        elif arg in ['dna', 'multiline']:
            result.append([x.strip() for x in sys.stdin.readlines()])

    if return_list:
        debug("Reading input:", result)
        return result
    else:
        debug("Reading input:", result[0])
        return result[0]

# End of GetInput()

def CompareResults(expected, found, ordered=True):

    if ordered:
        if expected == found:
            print("Results match")
        else:
            print("Results don't match")
            print("Expected:", expected)
            print("Found   :", found)
    else:
        if type(expected) != list:
            print("Results don't match.  Expected results are not a list.")
        elif type(found) != list:
            print("Results don't match.  Found results are not a list.")
        else:
            for exp in expected:
                if exp not in found:
                    print("Results don't match.", exp, "not found in result.")
                    return
            for fou in found:
                if fou not in expected:
                    print("Results don't match.", fou, "not expected.")
                    return
            print("Results match")

# End of CompareResult()

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

def FindPatternLocations(sequence, pattern, max_hamming_distance=0):

    locations = []
    len_sequence = len(sequence)
    len_pattern = len(pattern)
    srch_count = len_sequence - len_pattern + 1

    for i in range(0, srch_count):
        if FindHammingDistance(pattern, sequence[i:i+len_pattern]) <= max_hamming_distance:
            locations.append(i)

    return locations

# End of FindPatternLocations()

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
NumberToBaseComplementArray = ['T', 'G', 'C', 'A']

def NumberToPattern(number, k):

    result = ''
    for i in range(2*k-2, -2, -2):
        window = 0b11 << i
        result += NumberToBaseArray[(number & window)>>i]

    return result

# End of NumberToPattern()

def KmersFromText(text, k):

    result = []
    for i in range(len(text) - k + 1):
        result.append(text[i:i+k])

    return result

# End of KmersFromText()

def RandomWeighted(probabilities):

    cdf = [0] * len(probabilities)
    running_sum = 0
    for i, prob in enumerate(probabilities):
        running_sum += prob
        cdf[i] = running_sum

    rand = random.random() * running_sum

    # binary search
    upper_bound = len(cdf)-1
    lower_bound = 0
    while True:

        i = (upper_bound + lower_bound) // 2

        if rand < cdf[i]:
            upper_bound = i
        elif rand > cdf[i]:
            lower_bound = i + 1
        else:
            return i  # For the off-chance of an exact match

        if not (lower_bound < upper_bound):
            return upper_bound
    
# End of RandomWeighted()

def FindMostFrequentUsingDict(sequence, k, t=0):

    '''
    If t is non-zero, find all patterns with frequency greater than t.
    Otherwise, find the pattern that occurs the most.
    '''

    pattern_counts = OrderedDict()
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
    for pattern, count in pattern_counts.items():
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

def ComputeFrequencies(text, k):

    """Return a dictionary of k-mers to k-mer counts for the given text."""

    frequencies = {}

    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        if kmer in frequencies:
            frequencies[kmer] += 1
        else:
            frequencies[kmer] = 1

    return frequencies

# End of ComputeFrequencies()

def ReverseComplement(pattern):

    result = ''
    for ch in pattern:
        result += NumberToBaseComplementArray[BaseToNumber(ch)]
    result = result[::-1]   # Reverse the string
    return result

# End of ReverseComplement()

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
    window_end = L
    search_end = len(genome)
    patterns_found = []

    while window_end <= search_end:
        new_patterns = FindMostFrequentUsingDict(genome[window_start:window_end], k, t)
        for pattern in new_patterns:
            if pattern not in patterns_found:
                patterns_found.append(pattern)

        window_start += 1
        window_end += 1

    return patterns_found

# End of FindPatternClump()

def BetterClumpFinding(k, L, t, genome):

    clumped_patterns = []
    pattern_counts = ComputeFrequencies(genome[0:L], k)
    for pattern, count in pattern_counts.items():
        if count >= t:
            # At this point we can assume that pattern is not already
            # in clumped_patterns
            clumped_patterns.append(pattern)

    for i in range(1,len(genome)-L+1):
        first = genome[i-1:i-1+k]
        last  = genome[i+L-k:i+L]
        pattern_counts[first] -= 1    # first is guaranteed to already
                                         # be in pattern_counts
        if last in pattern_counts:
            pattern_counts[last] += 1
        else:
            pattern_counts[last] = 1
        if pattern_counts[last] >= t and last not in clumped_patterns:
            clumped_patterns.append(last)

    return sorted(clumped_patterns)

# End of BetterClumpFinding()

def Skew(genome):

    cur_skew = 0
    skew_diagram = [0]

    for i, base in enumerate(genome):
        if base == 'G':
            cur_skew += 1
        elif base == 'C':
            cur_skew -= 1

        skew_diagram.append(cur_skew)

    return skew_diagram

# End of Skew()

def FindMinimumSkews(genome):

    cur_skew = 0
    min_skew = 0
    min_skew_locci = []

    # debug_skews = []

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

        # debug_skews.append(cur_skew)

    return min_skew_locci

# End of FindMinimumSkews()

def FindHammingDistance(pattern_a, pattern_b):

    hamming_distance = 0

    for i in range(0, len(pattern_a)):
        if pattern_a[i] != pattern_b[i]:
            hamming_distance += 1

    return hamming_distance

# End of FindHammingDistance()

def DistanceBetweenPatternAndStrings(pattern, dna):

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

    return hamming_sum

# End of DistanceBetweenPatternAndStrings()

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

    patterns_found = []
    max_count_found = 0

    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        pattern_count = PatternCount(text, pattern, d)
        if pattern_count == max_count_found:
            patterns_found.append(pattern)
        elif pattern_count > max_count_found:
            patterns_found = [pattern]
            max_count_found = pattern_count

    return patterns_found

# End of FindMostFrequentWithMismatches()

def FindMostFrequentWithMismatchesAndReverseComplement(k, d, text):

    '''Look for all possible patterns in the given text'''

    patterns_found = []
    max_count_found = 0

    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        reverse_complement = ReverseComplement(pattern)
        pattern_count = PatternCount(text, pattern, d)
        pattern_count += PatternCount(text, reverse_complement, d)
        if pattern_count == max_count_found:
            patterns_found.append(pattern)
        elif pattern_count > max_count_found:
            patterns_found = [pattern]
            max_count_found = pattern_count

    return patterns_found

# End of FindMostFrequentWithMismatchesAndReverseComplement()

def FindMostFrequentWithMismatchesAndReverseComplementWithLocci(k, d, text):

    '''Look for all possible patterns in the given text'''

    patterns_found = {}
    max_count_found = 2   # Handle case where no patterns are found

    progress_divisor = 16
    progress_interval = int((4**k)/progress_divisor)
    progress_marker = 0

    for i in range(4**k):
        if (i%progress_interval) == 0:
            progress_marker += 1
        pattern = NumberToPattern(i, k)
        reverse_complement = ReverseComplement(pattern)
        pattern_locations = FindPatternLocations(text, pattern, d)
        pattern_locations += FindPatternLocations(text, reverse_complement, d)
        pattern_count = len(pattern_locations)
        if pattern_count == max_count_found:
            patterns_found[pattern] = pattern_locations
        elif pattern_count > max_count_found:
            patterns_found = {pattern : pattern_locations}
            max_count_found = pattern_count

    return patterns_found

# End of FindMostFrequentWithMismatchesAndReverseComplementWithInfo()



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

def MotifEnumeration(k, d, dna):

    pattern_lists = []

    # Find all of the patterns per dna strand
    for strand in dna:
        pattern_list = []
        for i in range(len(strand)-k+1):
            pattern = strand[i:i+k]
            neighbors = FindNeighbors(pattern, d)
            for neighbor in neighbors:
                if neighbor not in pattern_list:
                    pattern_list.append(neighbor)
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
        if pattern_in_all_lists:
            unique_patterns.append(pattern)

    return unique_patterns

# End of MotifEnumeration()

def MotifEnumeration2(k, d, dna):

    patterns = set()

    for strand in dna:
        for kmer in KmersFromText(strand, k):
            for neighbor in FindNeighbors(kmer, d):
                pattern_found_in_all_strands = True
                for strand2 in dna:
                    if PatternCount(strand2, neighbor, d) == 0:
                        pattern_found_in_all_strands = False
                        break
                if pattern_found_in_all_strands:
                    patterns.add(neighbor)

    return list(patterns)

# End of MotifEnumeration2()

def FindMedianString(k, dna):

    t = len(dna)
    distance = k * t
    median = ''
    saved_median_distances = OrderedDict()
    saved_medians_of_min_distance = []

    # Find distance for each possible k-mer
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        new_distance = DistanceBetweenPatternAndStrings(pattern, dna)
        #print ("New hamming sum is:", new_distance)
        saved_median_distances[pattern] = new_distance
        if new_distance < distance:
            distance = new_distance
            median = pattern
            saved_medians_of_min_distance = [pattern]
        elif new_distance == distance:
            saved_medians_of_min_distance.append(pattern)

    return median

# End of FindMedianString()

def GenProfileProbabilitiesKmer(kmer, profile):

    return reduce(lambda x, y: x*y,
               [profile[BaseToNumber(base)][i] for i, base in enumerate(kmer)])

# End of GenProfileProbabilitiesKmer()

def ProfileMostProbableKmer(text, k, profile):

    result = text[:k]   # Pick this if all k-mers have 0 probability
    best_probability = 0

    for i in range (len(text) - k + 1):
        kmer = text[i:i+k]
        cur_probability = GenProfileProbabilitiesKmer(kmer, profile)

        if cur_probability > best_probability:
            best_probability = cur_probability
            result = kmer

    return result

# End of ProfileMostProbableKmer()

def ProfileRandomlyGeneratedKmer(text, k, profile):

    result = text[:k]   # Pick this if all k-mers have 0 probability
    probabilities = []

    for i in range (len(text) - k + 1):
        kmer = text[i:i+k]
        cur_probability = GenProfileProbabilitiesKmer(kmer, profile)
        probabilities.append(cur_probability)

    result_i = RandomWeighted(probabilities)
    result = text[result_i:result_i+k]

    return result

# End of ProfileRandomlyGeneratedKmer()

def GenProfileFromMotifsBasic(motifs):

    num_motifs = len(motifs)
    motif_len = len(motifs[0])
    # profile = list(itertools.repeat(list(itertools.repeat(0.0, motif_len)), 4))
    profile = [
        list(itertools.repeat(0, motif_len)),
        list(itertools.repeat(0, motif_len)),
        list(itertools.repeat(0, motif_len)),
        list(itertools.repeat(0, motif_len))
    ]

    for i in range(num_motifs):
        for j in range(motif_len):
            profile[BaseToNumber(motifs[i][j])][j] += 1.0/num_motifs

    return profile

# End of GenProfileFromMotifsBasic()

def GenProfileFromMotifs(motifs):

    num_motifs = len(motifs)
    divisor = num_motifs + 4    # Add one for each base
    init_val = 1.0/divisor
    motif_len = len(motifs[0])
    profile = [
        [init_val] * motif_len,
        [init_val] * motif_len,
        [init_val] * motif_len,
        [init_val] * motif_len
    ]

    for i in range(num_motifs):
        for j in range(motif_len):
            profile[BaseToNumber(motifs[i][j])][j] += 1.0/divisor

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
        score += num_motifs - max(counts)

    return score

# End of CalculateScore()

def GreedyMotifSearch(k, t, dna):

    best_motifs = [x[:k] for x in dna]

    for j in range(len(dna[0]) - k + 1):
        motifs = [dna[0][j:j+k]]

        for i in range(1, t):

            profile = GenProfileFromMotifs(motifs)
            motifs.append(ProfileMostProbableKmer(dna[i], k, profile))

        if CalculateScore(motifs) < CalculateScore(best_motifs):
            best_motifs = motifs

    return best_motifs

# End of GreedyMotifSearch()

def Motifs(profile, dna):

    t = len(dna)
    k = len(profile[0])
    motifs = []
    for i in range(t):
        motifs.append(ProfileMostProbableKmer(dna[i], k, profile))

    return motifs

# End of Motifs()

def RandomizedMotifSearch(k, t, dna):

    # Generate initial random motifs
    len_string = len(dna[0])
    motifs = []
    for j in range(t):
        rand_num = random.randrange(len_string-k+1)
        motifs.append(dna[j][rand_num:rand_num+k])

    best_motifs = motifs

    while True:
        profile = GenProfileFromMotifs(motifs)
        motifs = Motifs(profile,dna)

        if CalculateScore(motifs) < CalculateScore(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs

# End of RandomizedMotifSearch()

def RandomizedMotifSearchBestOfMany(k, t, dna, iter_count):

    best_score = k*t
    best_result = None
    for i in range(iter_count):
        result = RandomizedMotifSearch(k, t, dna)
        cur_score = CalculateScore(result)
        if cur_score < best_score:
            debug(cur_score, "<--", ' '.join(result))
            best_score = cur_score
            best_result = result

    return best_result

# End of RandomizedMotifSearchBestOfMany()

def GibbsSampler(k, t, N, dna):

    # Generate initial random motifs
    len_string = len(dna[0])
    motifs = []
    for j in range(t):
        rand_num = random.randrange(len_string-k+1)
        motifs.append(dna[j][rand_num:rand_num+k])

    best_motifs = motifs
    best_score = t * k

    for j in range(N):

        i = random.randrange(t)
        limited_motifs = motifs[:i] + motifs[i+1:]
        profile = GenProfileFromMotifs(limited_motifs)
        new_motif = ProfileRandomlyGeneratedKmer(dna[i], k, profile)
        motifs = motifs[:i] + [new_motif] + motifs[i+1:]

        cur_score = CalculateScore(motifs)
        if cur_score < best_score:
            debug(cur_score, "<--", ' '.join(motifs))
            best_motifs = motifs
            best_score = cur_score

    return best_motifs

# End of GibbsSampler()

def Composition(k, text):

    result = []

    for i in range(len(text) - k + 1):
        result.append(text[i:i+k])

    return sorted(result)

# End of Composition()

def StringFormedByGenomePath(genome_path, start_prefix = None):

    k = len(genome_path[0])
    N = len(genome_path)
    new_string = genome_path[0]
    for i in range(1,N):
        new_string += genome_path[i][k-1]

    return new_string

# End of StringFormedByGenomePath()

def Overlap(patterns):

    result = []  # A list of lists.  The first element of an inner list contains
                 # a pattern, and all subsequent elements of an inner list
                 # are potential downstream neighbors of the first element.
    for i, pat in enumerate(patterns):
        suffix = pat[1:]
        remaining_patterns = list(patterns)
        remaining_patterns.remove(pat)
        result.append([pat])
        for rpat in remaining_patterns:
            if rpat.startswith(suffix):
                result[i].append(rpat)

    return result

# End of Overlap()

def StitchPatterns(patterns, cur_path = [], solutions = [], start_prefix = None):

    if len(patterns) == 0:
        solutions.append(cur_path)

    pool_of_next = []
    if start_prefix != None:
        # We're restricting the possible next option based on a starting prefix.
        for pattern in patterns:
            if pattern.startswith(start_prefix):
                pool_of_next.append(pattern)
    else:
        # Use the entire genome path
        pool_of_next = list(patterns)

    if len(pool_of_next) == 0:
        # Abandon this path
        return

    for next in pool_of_next:
        new_start_prefix = next[1:]
        new_patterns = list(patterns)
        new_patterns.remove(next)
        StitchPatterns(new_patterns, cur_path + [next], solutions, new_start_prefix)

    return solutions

# End of StitchPatterns()

class Node:

    def __init__(self, label):

        self.label = label
        self.outgoing_edges = []
        self.incoming_edges = []
        self.indegree = 0
        self.outdegree = 0

    # end of Node.__init__()

    def __str__(self):
        return str(self.label)
    # End of Node.__str__()

    def next_node(self):
        # Assumes a linear sequence of nodes
        if len(self.outgoing_edges) > 0:
            return self.outgoing_edges[0].next_node
        else:
            return None
    # End of Node.next_node()

    def next_nodes(self):
        if len(self.outgoing_edges) > 0:
            return sorted([x.next_node for x in self.outgoing_edges])
        else:
            return []
    # End of Node.next_nodes()

    def find_unvisited_edge(self):
        for edge in self.outgoing_edges:
            if not edge.visited:
                return edge
        return None
    # End of Node.find_unused_edge()

    def __lt__(self, other):
        return self.label < other.label
    # End of Node.__lt__()

    def __eq__(self, other):
        return self.label == other.label
    # End of Node.__eq__()

# End of class Node

class Edge:

    def __init__(self, label=None, edge=None):
        if label is not None:
            self.label = label
        elif edge is not None:
            self.label = edge.label
        else:
            raise Exception("Missing label or edge parameter")
        self.prev_node = None
        self.next_node = None
        self.visited = False
    # End of Edge.__init__()

    def __str__(self):
        return self.label
    # End of Edge.__str__()

# End of class Edge

class PathGraphk:

    def __init__(self, k, text):

        self.head_node = None

        edges = [Edge(text[i:i+k]) for i in range(len(text)-k+1)]
        prev_edge = None
        for edge in edges:
            # Create the node before the current edge
            node = Node(edge.label[:-1])
            node.outgoing_edges.append(edge)
            if self.head_node is None:
                # Remember the head node
                self.head_node = node
            if prev_edge:
                prev_edge.next_node = node
            prev_edge = edge

        # Add the last node
        node = Node(prev_edge.label[1:])
        prev_edge.next_node = node

    # End of PathGraphk.__init__()

    def __iter__(self):
        self.cur_node = self.head_node
        return self
    # End of PathGraphk.__iter__()

    def __next__(self):
        ret_node = self.cur_node
        if ret_node is None:
            raise StopIteration
        self.cur_node = self.cur_node.next_node()
        return ret_node
    # End of PathGraphk.__next__()

# End of class PathGraphk

class DeBruijn_Parent(object):

    def __init__(self):

        self.node_map = {}
        self.node_list = []
        self.edge_list = []
        self.all_graphs = []
        self.k = 0
        self.d = 0
        self.kd_graph = False
        self.e_cycle = []

    # End of DeBruijn_Parent.__init__()

    def add_node(self, label):
        hash_key = str(label)
        if hash_key in self.node_map:
            debug("Merging label:", label)
            node = self.node_map[hash_key]
        else:
            debug("Adding  label:", label)
            node = Node(label)
            self.node_map[hash_key] = node
            self.node_list.append(node)
        return node
    # End of DeBruijn_Parent.add_node()

    def EulerianCycle(self):

        if len(self.node_map) == 0:
            return None

        self.reset_cycles()
        cur_node = list(self.node_map.values())[0]  # Pick a start node at random
        debug("Starting at node: ", cur_node)
        cur_cycle = [cur_node]

        while True:   # while unvisited edges remain

            # For each edge in graph
            while True:
                next_edge = cur_node.find_unvisited_edge()
                if next_edge is not None:
                    next_edge.visited = True
                    cur_node = next_edge.next_node
                    cur_cycle.append(cur_node)
                else:
                    break

            #time.sleep(2)
            debug('->'.join([str(x) for x in cur_cycle]))
            cycle_index = self.find_cycle_index_of_node_with_unused_edges(cur_cycle)

            if cycle_index is -1:
                # Every edge has been visited, return last cycle found
                return cur_cycle
            else:
                # There are still unused edges in the graph, so make
                # the edge with unused edges the new starting point,
                # then keep going
                prev_cycle = cur_cycle
                cur_cycle = prev_cycle[cycle_index:]
                cur_cycle += prev_cycle[1:cycle_index]
                cur_node = cur_cycle[0]
                cur_cycle.append(cur_node)  
                debug("Restarting from node: ", cur_node)

        # End of while unvisited edges remain

        # We should never get here, but just in case...
        raise Exception("Something went wrong :(")

    # End of DeBruijn_Parent.EulerianCycle()

    def EulerianPath(self):

        # Find start and end nodes
        start_node = None
        end_node = None
        for node in self:
            if node.indegree == node.outdegree - 1:
                if start_node is None:
                    start_node = node
                else:
                    raise Exception("Multiple start nodes found: " +
                                    str(start_node) + " and " + str(node))
            elif node.outdegree == node.indegree -1:
                if end_node is None:
                    end_node = node
                else:
                    raise Exception("Multiple end nodes found: " +
                                    str(end_node) + " and " + str(node))
            elif node.outdegree != node.indegree:
                raise Exception("Invalidly unbalanced node: " + str(node) +
                                ".  Indegree = " + str(indegree) +
                                ".  Outdegree = " + str(outdegree) + ".")

        # Draw an edge from the end node to the start node to balance the graph
        if start_node is not None and end_node is not None:
            new_edge = Edge('to ' + str(start_node))
            end_node.outgoing_edges.append(new_edge)
            new_edge.next_node = start_node
            start_node.indegree += 1
            end_node.outdegree += 1
        elif start_node is None and end_node is None:
            # This is okay, it just means the start and end nodes were merged.
            pass
        else:
            raise Exception("Invalid number of unbalanced nodes for an Eulerian path.")

        # Now find an Eulerian cycle
        cycle = self.EulerianCycle()

        # Almost done.  Now we just have to shift the cycle to start at our
        # start node, and break the cycle into a path.
        if start_node is not None:
            for i in range(len(cycle)-1):
                j = i + 1
                if cycle[i] == end_node and cycle[j] == start_node:
                    new_cycle = cycle[j:-1]   # Chop off the loop back to start
                    new_cycle += cycle[:i+1]  # Capture cycle[i] here
                    cycle = new_cycle
                    break
        else:
            # If no recogniseable start/end nodes were found, just break the
            # cycle to make a path at its current end
            cycle.pop()

        if start_node is not None and cycle[0] != start_node:
            raise Exception("Start node of '" + str(start_node) + "' doesn't match actual start of cycle '"
                            + str(cycle[0]) + ".  Full cycle is: " + str([str(x) for x in cycle]))

        return cycle

    # End of DeBruijn_Parent.EulerianPath()

    def ReconstructPathString(self):

        path = self.EulerianPath()
        debug("Path is:", '->'.join([str(x) for x in path]))

        # Now stitch together the nodes
        new_string = ''
        for node in path:
            debug("Adding", node.label[0], "from", node.label, "to", new_string)
            new_string += node.label[0]
        debug("Finally adding", node.label[1:], "from", node.label, "to", new_string)
        new_string += node.label[1:]

        return new_string

    # End of DeBruijn_Parent.ReconstructPathString()

    def ReconstructCycleString(self):

        path = self.EulerianCycle()
        debug("Path is:", '->'.join([str(x) for x in path]))

        # Now stitch together the nodes
        new_string = ''
        for node in path[:-1]:
            # The last node is the same as the first (it's a cycle), so don't
            # take that one.
            debug("Adding", node.label[0], "from", node.label, "to", new_string)
            new_string += node.label[0]

        return new_string

    # End of DeBruijn_Parent.ReconstructCycleString()

    def FindMaximalNonBranchingPaths(self):

        paths = []

        for node in self.node_list:
            if ((node.indegree != 1 or node.outdegree != 1)
                and node.outdegree > 0):
                for edge in node.outgoing_edges:
                    non_branching_path = [edge]
                    edge.visited = True
                    next_node = edge.next_node
                    while next_node.indegree == 1 and next_node.outdegree == 1:
                        next_edge = next_node.outgoing_edges[0]
                        non_branching_path.append(next_edge)
                        next_edge.visited = True
                        next_node = next_edge.next_node
                    paths.append(non_branching_path)

        # Find isolated cycles, and add to paths
        for edge in self.edge_list:
            if not edge.visited:
                # All remaining unvisited edges will be for 1-in-1-out nodes.
                non_branching_path = [edge]
                edge.visited = True
                next_edge = edge.next_node.outgoing_edges[0]
                while not next_edge.visited:
                    non_branching_path.append(next_edge)
                    next_edge.visited = True
                    next_edge = next_edge.next_node.outgoing_edges[0]
                paths.append(non_branching_path)

        return paths

    # End of DeBruijn_Parent.FindMaximalNonBranchingPaths()

    def GetStringSpelledByPatterns(self, patterns=None):

        if patterns is None:
            patterns = [x.label for x in self.node_list]

        result = ''
        for pattern in patterns:
            result += pattern[0]
        result += pattern[1:]

        return result

    # End of DeBruijn_Parent.GetStringSpelledByPatterns()

    def GenContigs(self):

        max_non_branching_paths = self.FindMaximalNonBranchingPaths()
        contigs = [self.GetStringSpelledByPatterns([str(x) for x in path]) for path in max_non_branching_paths]
        return contigs

    # End of DeBruijn_Parent.GenContigs()

    def reset_cycles(self):
        for label in self.node_map:
            node = self.node_map[label]
            for edge in node.outgoing_edges:
                edge.visited = False
    # End of DeBruijn_Parent.reset_cycles()

    def find_next_node(self):

        cur_node = self.e_cycle[-1]
        next_node = None
        lead_node_index = len(self.e_cycle)
        minus_d_index = lead_node_index - self.d
        if self.kd_graph and minus_d_index >= 0:
            minus_d_node = self.e_cycle[minus_d_index]
        else:
            minus_d_node = None
        for edge in cur_node.outgoing_edges:
            if not edge.visited:
                found = False
                try_node = edge.next_node
                if minus_d_node is not None:
                    if try_node.label[0][-1] == minus_d_node.label[1][0]:
                        found = True
                else:
                    found = True
                if found:
                    edge.visited = True
                    next_node = try_node
                    break

        return next_node

    # End of DeBruijn_Parent.find_next_node()

    def find_cycle_index_of_node_with_unused_edges(self, cycle):
        for i, node in enumerate(cycle):
            unvisited_edge = node.find_unvisited_edge()
            if unvisited_edge is not None:
                return i
        return -1
    # End of DeBruijn_Parent.find_node_with_unused_edges()

    def __iter__(self):
        self.node_list = sorted(list(self.node_map.values()))
        self.i = 0
        return self
    # End of DeBruijn_Parent.__iter__()

    def __next__(self):
        if self.i >= len(self.node_list):
            raise StopIteration
        return_node = self.node_list[self.i]
        self.i += 1
        return return_node
    # End of DeBruijn_Parent.__next__()

# End of class DeBruijn_Parent

class DeBruijnk(DeBruijn_Parent):

    def __init__(self, *args):

        super().__init__()

        debug("In DeBruijnk.__init__().  Args are:", args)

        if (len(args) == 2 and
            type(args[0]) == int and
            type(args[1]) == str):
            self.__init_from_k_text(*args)
        elif (len(args) == 1 and
              type(args[0]) == list and
              type(args[0][0]) == str and
              ' -> ' in args[0][0]):
            self.__init_from_adjacency_list(*args)
        elif (len(args) == 1 and
              type(args[0]) == list and
              type(args[0][0]) == str):
            self.__init_from_patterns(*args)
        else:
            raise Exception("Invalid arguments: " + str(args))

    # End of DeBruijnk.__init__()

    def __init_from_k_text(self, k, text):

        edges = [Edge(text[i:i+k]) for i in range(len(text)-k+1)]
        debug("Edges:", [str(x) for x in edges])
        prev_edge = None
        for edge in edges:
            # Create the node before the current edge
            node = self.add_node(edge.label[:-1])
            node.outgoing_edges.append(edge)
            if prev_edge is not None:
                prev_edge.next_node = node
            prev_edge = edge

        # Add the last node
        if prev_edge is not None:
            node = self.add_node(edge.label[1:])
            prev_edge.next_node = node

    # End of DeBruijnk.__init_from_k_text()

    def __init_from_patterns(self, patterns):

        debug("Patterns are:", patterns)
        edges = [Edge(x) for x in patterns]

        for edge in edges:
            start_node = self.add_node(edge.label[:-1])
            end_node = self.add_node(edge.label[1:])
            start_node.outgoing_edges.append(edge)
            edge.next_node = end_node
            start_node.outdegree += 1
            end_node.indegree += 1

    # End of DeBruijnk.__init_from_patterns()

    def __init_from_adjacency_list(self, adjacency_list):

        for adjacency in adjacency_list:
            start_name, end_names = adjacency.split(' -> ')
            end_names = end_names.split(',')
            start_node = self.add_node(start_name)
            for end_name in end_names:
                end_node = self.add_node(end_name)
                new_edge = Edge(label=start_name+' to '+end_name)
                self.edge_list.append(new_edge)
                start_node.outgoing_edges.append(new_edge)
                end_node.incoming_edges.append(new_edge)
                new_edge.next_node = end_node
                new_edge.prev_node = start_node
                start_node.outdegree += 1
                end_node.indegree += 1

    # End of DeBruijnk.__init_from_adjacency_list()

#------------------------------------------------------------------------------------
#    def EulerianCycle(self):                                                        
#                                                                                    
#        if len(self.node_map) == 0:                                                 
#            return None                                                             
#                                                                                    
#        self.reset_cycles()                                                         
#        cur_node = list(self.node_map.values())[0]  # Pick a start node at random   
#        debug("Starting at node: ", cur_node)                                       
#        cur_cycle = [cur_node]                                                      
#                                                                                    
#        while True:   # while unvisited edges remain                                
#                                                                                    
#            # For each edge in graph                                                
#            while True:                                                             
#                next_edge = cur_node.find_unvisited_edge()                          
#                if next_edge is not None:                                           
#                    next_edge.visited = True                                        
#                    cur_node = next_edge.next_node                                  
#                    cur_cycle.append(cur_node)                                      
#                else:                                                               
#                    break                                                           
#                                                                                    
#            #time.sleep(2)                                                          
#            debug('->'.join([str(x) for x in cur_cycle]))                           
#            cycle_index = self.find_cycle_index_of_node_with_unused_edges(cur_cycle)
#                                                                                    
#            if cycle_index is -1:                                                   
#                # Every edge has been visited, return last cycle found              
#                return cur_cycle                                                    
#            else:                                                                   
#                # There are still unused edges in the graph, so make                
#                # the edge with unused edges the new starting point,                
#                # then keep going                                                   
#                prev_cycle = cur_cycle                                              
#                cur_cycle = prev_cycle[cycle_index:]                                
#                cur_cycle += prev_cycle[1:cycle_index]                              
#                cur_node = cur_cycle[0]                                             
#                cur_cycle.append(cur_node)                                          
#                debug("Restarting from node: ", cur_node)                           
#                                                                                    
#        # End of while unvisited edges remain                                       
#                                                                                    
#        # We should never get here, but just in case...                             
#        raise Exception("Something went wrong :(")                                  
#                                                                                    
#    # End of DeBruijnk.EulerianCycle()                                              
#------------------------------------------------------------------------------------

    def EulerianPath(self):

        # Find start and end nodes
        start_node = None
        end_node = None
        for node in self:
            if node.indegree == node.outdegree - 1:
                if start_node is None:
                    start_node = node
                else:
                    raise Exception("Multiple start nodes found: " +
                                    str(start_node) + " and " + str(node))
            elif node.outdegree == node.indegree -1:
                if end_node is None:
                    end_node = node
                else:
                    raise Exception("Multiple end nodes found: " +
                                    str(end_node) + " and " + str(node))
            elif node.outdegree != node.indegree:
                raise Exception("Invalidly unbalanced node: " + str(node) +
                                ".  Indegree = " + str(indegree) +
                                ".  Outdegree = " + str(outdegree) + ".")

        # Draw an edge from the end node to the start node to balance the graph
        if start_node is not None and end_node is not None:
            new_edge = Edge('to ' + str(start_node))
            end_node.outgoing_edges.append(new_edge)
            new_edge.next_node = start_node
            start_node.indegree += 1
            end_node.outdegree += 1
        elif start_node is None and end_node is None:
            # This is okay, it just means the start and end nodes were merged.
            pass
        else:
            raise Exception("Invalid number of unbalanced nodes for an Eulerian path.")

        # Now find an Eulerian cycle
        cycle = self.EulerianCycle()

        # Almost done.  Now we just have to shift the cycle to start at our
        # start node, and break the cycle into a path.
        if start_node is not None:
            for i in range(len(cycle)-1):
                j = i + 1
                if cycle[i] == end_node and cycle[j] == start_node:
                    new_cycle = cycle[j:-1]   # Chop off the loop back to start
                    new_cycle += cycle[:i+1]  # Capture cycle[i] here
                    cycle = new_cycle
                    break
        else:
            # If no recogniseable start/end nodes were found, just break the
            # cycle to make a path at its current end
            cycle.pop()

        if start_node is not None and cycle[0] != start_node:
            raise Exception("Start node of '" + str(start_node) + "' doesn't match actual start of cycle '"
                            + str(cycle[0]) + ".  Full cycle is: " + str([str(x) for x in cycle]))

        return cycle

    # End of DeBruijnk.EulerianPath()

    def ReconstructPathString(self):

        path = self.EulerianPath()
        debug("Path is:", '->'.join([str(x) for x in path]))

        # Now stitch together the nodes
        new_string = ''
        for node in path:
            debug("Adding", node.label[0], "from", node.label, "to", new_string)
            new_string += node.label[0]
        debug("Finally adding", node.label[1:], "from", node.label, "to", new_string)
        new_string += node.label[1:]
        
        return new_string

    # End of DeBruijnk.ReconstructPathString()

    def ReconstructCycleString(self):

        path = self.EulerianCycle()
        debug("Path is:", '->'.join([str(x) for x in path]))

        # Now stitch together the nodes
        new_string = ''
        for node in path[:-1]:
            # The last node is the same as the first (it's a cycle), so don't
            # take that one.
            debug("Adding", node.label[0], "from", node.label, "to", new_string)
            new_string += node.label[0]

        return new_string

    # End of DeBruijnk.ReconstructCycleString()

# End of class DeBruijnk

class PairedDeBruijn(DeBruijn_Parent):

    def __init__(self, k, d, kdmers):

        super().__init__()

        self.k = k
        self.d = d
        self.kdmers = kdmers
        self.kd_graph = True

        debug("(" + str(k) + "," + str(d) + ")-mers are:", kdmers)
        self.edge_list = [Edge(x) for x in kdmers]

        for edge in self.edge_list:
            p1,p2 = str(edge).split('|')
            start_label = [p1[:-1], p2[:-1]]
            end_label   = [p1[1:],  p2[1:]]
            start_node = self.add_node(start_label)
            end_node = self.add_node(end_label)
            start_node.outgoing_edges.append(edge)
            end_node.incoming_edges.append(edge)
            edge.prev_node = start_node
            edge.next_node = end_node
            start_node.outdegree += 1
            end_node.indegree += 1

    # End of PairedDeBruijn.__init__()

    def ReconstructPathString(self):

        path = super().EulerianPath()
        debug("Path is:", '->'.join([str(x) for x in path]))

        # Now stitch together the nodes

        first_patterns = []
        second_patterns = []
        for node in path:
            first, second = node.label
            first_patterns.append(first)
            second_patterns.append(second)

        prefix_string = self.GetStringSpelledByPatterns(first_patterns)
        suffix_string = self.GetStringSpelledByPatterns(second_patterns)

        for i in range(self.k + self.d, len(prefix_string)):
            if prefix_string[i] != suffix_string[i-(self.k+self.d)]:
                return "there is no string spelled by the gapped patterns"

        result = prefix_string + suffix_string[-(self.k+self.d):]

        return result


#        new_string = ''
#        for node in path:
#            debug("Adding", node.label[0][0], "from", node.label, "to", new_string)
#            new_string += node.label[0][0]
#        debug("Finally adding", node.label[1][1:], "from", node.label, "to", new_string)
#        new_string += node.label[1][1:]
#
#        return new_string

    # End of PairedDeBruijn.ReconstructPathString()

    def GetStringSpelledByGappedPatterns(self):

        first_patterns = []
        second_patterns = []
        for edge in self.edge_list:
            first, second = edge.label.split('|')
            first_patterns.append(first)
            second_patterns.append(second)

        prefix_string = self.GetStringSpelledByPatterns(first_patterns)
        suffix_string = self.GetStringSpelledByPatterns(second_patterns)

        for i in range(self.k + self.d, len(prefix_string)):
            if prefix_string[i] != suffix_string[i-(self.k+self.d)]:
                return "there is no string spelled by the gapped patterns"

        result = prefix_string + suffix_string[-(self.k+self.d):]

        return result

    # End of PairedDeBruijn.GetStringSpelledByGappedPatterns()

# End of class PairedDeBruijn


def KUniversalCircularString(k):

    patterns = ["{x:0{k}b}".format(k=k,x=x) for x in range(2**k)]

    graph = DeBruijnk(patterns)

    new_string = graph.ReconstructCycleString()

    return new_string

# End of KUniversalCircularString()

class AllGraphs(object):

    def __init__(self, start_graph):

        self.all_graphs = [start_graph]

    # End of AllGraphs.__init__()

    def SplitIntoSimpleDirectedGraphs(self):

        while self.__find_complex_node():
            # Values from __find_complex_node are returned through member variables
            graph = self.complex_graph
            node_index  = self.complex_node_index
            node = graph.node_list[node_index]

            for incoming_edge_index in range(len(node.incoming_edges)):
                for outgoing_edge_index in range(len(node.outgoing_edges)):
                    new_graph = copy.deepcopy(graph)
                    node_to_split = new_graph.node_list[node_index]
                    incoming_edge_to_save = node_to_split.incoming_edges[incoming_edge_index]
                    outgoing_edge_to_save = node_to_split.outgoing_edges[outgoing_edge_index]
                    # remove all other incoming and outgoing edges
                    for incoming_edge in node_to_split.incoming_edges:
                        if incoming_edge is not incoming_edge_to_save:
                            self.__remove_edge(incoming_edge)
                    for outgoing_edge in node_to_split.outgoing_edges:
                        if outgoing_edge is not outgoing_edge_to_save:
                            self.__remove_edge(outgoing_edge)
                    self.all_graphs.append(new_graph)
            self.all_graphs.remove(graph)

    # End of SplitIntoSimpleDirected()

    def EliminateInvalidPairedPaths(self):

        pass

    # End of EliminateInvalidPairedPaths()

    def __find_complex_node(self):

        # Return data through class members, to allow method to
        # be called in a conditional.
        self.complex_graph = None
        self.complex_node_index = -1

        for graph in self.all_graphs:
            for i, node in enumerate(graph.node_list):
                if node.indegree > 1 or node.outdegree > 1:
                    self.complex_graph = graph
                    self.complex_node_index = i
                    return True

        return False

    # End of AllGraphs.__find_complex_node()

    def __remove_edge(self, edge):

        edge.prev_node.outgoing_edges.remove(edge)
        edge.prev_node.outdegree -= 1
        edge.next_node.incoming_edges.remove(edge)
        edge.next_node.indegree -= 1

    # End of __remove_edge()
                    
# End of AllGraphs()

RNACodonToAA = {'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAU':'N',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
                'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGU':'S',
                'AUA':'I', 'AUC':'I', 'AUG':'M', 'AUU':'I',

                'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAU':'H',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
                'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',

                'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAU':'D',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
                'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',

                'UAA':'*', 'UAC':'Y', 'UAG':'*', 'UAU':'Y',
                'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
                'UGA':'*', 'UGC':'C', 'UGG':'W', 'UGU':'C',
                'UUA':'L', 'UUC':'F', 'UUG':'L', 'UUU':'F'
                }

AAToRNACodons = {}
for codon, aa in RNACodonToAA.items():
    if aa in AAToRNACodons:
        AAToRNACodons[aa].append(codon)
    else:
        AAToRNACodons[aa] = [codon]

def TranslateRNAToProtein(rna):

    result = ''
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        aa = RNACodonToAA[codon]
        if aa != '*':
            result += aa
        else:
            break

    return result

# End of TranslateRNAToProtein()

def FindGeneForPeptide(text, peptide):

    results = find_gene_for_peptide_on_one_strand(text, peptide)

    rc_text = ReverseComplement(text)
    rc_results = find_gene_for_peptide_on_one_strand(rc_text, peptide)
    for res in rc_results:
        results.append(ReverseComplement(res))

    return results

# End of FindGeneForPeptide()

def find_gene_for_peptide_on_one_strand(text, peptide):

    results = []

    len_gene = len(peptide) * 3
    len_peptide = len(peptide)
    for i in range(0, len(text) - len_gene + 1):
        gene_candidate = text[i:i+len_gene]
        if does_gene_match_peptide(gene_candidate, peptide):
            results.append(gene_candidate)

    return results

# End of find_gene_for_peptide_in_one_reading_frame()

def does_gene_match_peptide(gene_candidate, peptide):

    result = True

    for peptide_index in range(len(peptide)):
        gene_index = peptide_index * 3
        test_codon = gene_candidate[gene_index:gene_index+3]
        test_rna_codon = test_codon.replace('T', 'U')
        test_aa = peptide[peptide_index]
        if len(test_rna_codon) != 3:
            print("Test RNA codon is not the proper length:", test_rna_codon)
            import pdb; pdb.set_trace()
        if RNACodonToAA[test_rna_codon] != test_aa:
            result = False
            break

    return result

# End of does_gene_match_peptide()

AAtoDaltons = {'G':  57, 'A':  71, 'S':  87, 'P':  97, 'V':  99, 'T': 101,
               'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
               'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
               'Y': 163, 'W': 186}

def GenTheoreticalSpectrum(peptide, cyclic=True):

    prefix_mass = [0]
    len_peptide = len(peptide)

    for i in range(len_peptide):
        aa = peptide[i]
        if type(aa) == str:
            aa_mass = AAtoDaltons[aa]
        else:
            aa_mass = aa
        prefix_mass.append(prefix_mass[-1] + aa_mass)

    cyclic_spectrum = [0]
    for i in range(len_peptide+1):
        for j in range(i+1, len_peptide+1):
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if cyclic and i > 0 and j < len_peptide:
                cyclic_spectrum.append(prefix_mass[len_peptide] - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spectrum)

# End of GenTheoreticalSpectrum()

def CountingPeptidesWithGivenMass(mass):

    #possible_masses = sorted(list(set(AAtoDaltons.values())))
    possible_masses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115,
                       128, 129, 131, 137, 147, 156, 163, 186]
    matches = []

    def recurse(remaining_mass, level = 0, aa_list = []):

        nonlocal matches

        debug("CountingPeptidesWithGivenMass.recurse():", level*'  ', "Mass:",
              remaining_mass, " Level:", level)

        count = 0
        for candidate_mass in possible_masses:
            #debug("CountingPeptidesWithGivenMass.recurse():", level*'  ',
            #      "Considering adding mass:", candidate_mass,
            #      "Remaining mass will be:", remaining_mass - candidate_mass)
            new_aa_list = aa_list + [candidate_mass]

            if candidate_mass == remaining_mass:
                debug("CountingPeptidesWithGivenMass.recurse():", level*'  ',
                      "Found match with:", new_aa_list)
                matches.append(new_aa_list)
                count += 1
                break
            elif candidate_mass < remaining_mass:
                new_remaining_mass = remaining_mass - candidate_mass
                sub_count = recurse(new_remaining_mass, level+1, new_aa_list)
                count += sub_count
            else:
                break # No more candidates are worth exploring

        return count

    # End of CountingPeptidesWithGivenMass.recurse()

    result = recurse(mass)

    debug("Found the following matches:\n", matches)

    return result

# End of CountingPeptidesWithGivenMass()

def ConsistentSpectrums(peptide_spectrum, gene_spectrum):
    '''Return true if the peptide spectrum is a subset of the gene spectrum'''

    # Make a copies of the spectra so that we don't destroy the originals
    peptide_spectrum = peptide_spectrum[:]
    gene_spectrum = gene_spectrum[:]

    for mass in peptide_spectrum:
        if mass in gene_spectrum:
            gene_spectrum.remove(mass)
        else:
            return False

    return True

# End of ConsistentSpectrums()

def CyclopeptideSequencing(spectrum):

    '''Each peptide is represented as a list of masses'''

    possible_masses = sorted(list(set(AAtoDaltons.values())))
    parent_mass = spectrum[-1]
    
    def expand(peptides):

        new_peptides = []
        for peptide in peptides:
            for mass in possible_masses:
                new_peptides.append(peptide + [mass])
        return new_peptides

    # End of CyclopeptideSequencing.expand()

    result = []
    peptides = [[]]

    while len(peptides) > 0:
        peptides = expand(peptides)
        debug("Expanded peptides to length", len(peptides[-1]))
        debug("Total # of peptides in current list:", len(peptides))
        for peptide in peptides[:]:
            if sum(peptide) == parent_mass:
                if GenTheoreticalSpectrum(peptide, cyclic=True) == spectrum:
                    debug("Adding:", peptide)
                    result.append(peptide)
                peptides.remove(peptide)
            elif not ConsistentSpectrums(
                     GenTheoreticalSpectrum(peptide, cyclic=False), spectrum):
                peptides.remove(peptide)
                    
    return result

# End of CyclopeptideSequencing()

def CyclopeptideScoring(peptide, spectrum):

    calculated_spectrum = GenTheoreticalSpectrum(peptide, cyclic=True)
    score = 0
    for mass in spectrum:
        if mass in calculated_spectrum:
            score += 1
            calculated_spectrum.remove(mass)  # Don't count it more than once

    return score

# End of CyclopeptideScoring()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Perform a bioinformatics assignment.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Turn on debug messages')
    parser.add_argument('-c', '--compare', action='store_true',
                        help='Compare results with expected results')
    args = parser.parse_args()

    if args.debug:
        DEBUG = True

    if True:

        #################################################################
        # Get input for Generate Contigs problem
        #################################################################

        kmers = GetInput('multiline')
        graph = DeBruijnk(kmers)
        result = graph.GenContigs()
        print(' '.join(result))

        if args.compare:
            expected = []
            with open('expected_output.txt', 'r') as er:
                for line in er.readlines():
                    expected += line.strip().split()
            CompareResults(expected, result, ordered=False)

    elif False:

        #################################################################
        # Get input for Linear Spectrum problem
        #################################################################

        peptide = GetInput('s')
        result = GenTheoreticalSpectrum(peptide, cyclic=False)
        print(' '.join([str(x) for x in result]))

        if args.compare:
            expected = []
            with open('expected_output.txt', 'r') as er:
                for line in er.readlines():
                    expected += [int(x) for x in line.strip().split()]
            CompareResults(expected, result, ordered=True)

    if False:

        #################################################################
        # Input for CyclopeptideSequencing problem
        #################################################################

        spectrum = GetInput('i i')
        result = CyclopeptideSequencing(spectrum)
        string_result = ' '.join(['-'.join([str(x) for x in peptide]) for peptide in result])
        print(string_result)

        if args.compare:
            expected = []
            with open('expected_output.txt', 'r') as er:
                for line in er.readlines():
                    expected += [[int(aa) for aa in peptide.split('-')] for peptide in line.strip().split()]
            CompareResults(expected, result, ordered=False)

    elif False:

        #################################################################
        # Saving the below to go back to the paired DeBruijn path problem
        #################################################################

        # DeBruijnk String Reconstruction Problem
        k, patterns = GetInput('i', 'multiline')
        start_time = time.time()
        graph = DeBruijnk(patterns)
        text_result = graph.ReconstructPathString()
        end_time = time.time()
        debug("Elapsed time is", end_time - start_time)

        string_result = text_result
        print(string_result)

        if args.compare:
            with open('expected_output.txt', 'r') as er:
                expected = [x.strip() for x in er.readlines()]
            CompareResults(expected[0], string_result, ordered=True)

    elif False:
       # Paired DeBruijnk String Reconstruction Problem
       k, d, patterns = GetInput('i i', 'multiline')
       start_time = time.time()
       graph = PairedDeBruijn(k, d, patterns)
       path_string = graph.ReconstructPathString()
       end_time = time.time()
       debug("Elapsed time is", end_time - start_time)

       string_result = path_string
       print(path_string)

       if args.compare:
           with open('expected_output.txt', 'r') as er:
               expected = [x.strip() for x in er.readlines()]
           CompareResults(expected[0], string_result, ordered=True)

    else:
        pass

