#!/usr/bin/env python

import sys
import array
import itertools
import random


def PatternCount(sequence, pattern, max_hamming_distance=0):

    count = 0
    len_sequence = len(sequence)
    len_pattern = len(pattern)
    srch_count = len_sequence - len_pattern + 1

    for i in range(0, srch_count):
        if FindHammingDistance(pattern, sequence[i:i+len_pattern]) <= max_hamming_distance:
            print "found match at i =", i
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

def KmersFromText(text, k):

    result = []
    for i in range(len(text) - k + 1):
        result.append(text[i:i+k])

    return result

# End of KmersFromText()

def RandomWeighted(probabilities):

    cdf = list(itertools.repeat(0,len(probabilities)))
    running_sum = 0
    for i, prob in enumerate(probabilities):
        running_sum += prob
        cdf[i] = running_sum

    rand = random.random() * running_sum

    # binary search
    upper_bound = len(cdf)-1
    lower_bound = 0
    while True:

        i = (upper_bound + lower_bound) / 2

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

def GenProfileFromMotifs(motifs):

    num_motifs = len(motifs)
    divisor = num_motifs + 4    # Add one for each base
    init_val = 1.0/divisor
    motif_len = len(motifs[0])
    # profile = list(itertools.repeat(list(itertools.repeat(0.0, motif_len)), 4))
    profile = [
        list(itertools.repeat(init_val, motif_len)),
        list(itertools.repeat(init_val, motif_len)),
        list(itertools.repeat(init_val, motif_len)),
        list(itertools.repeat(init_val, motif_len))
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

def GibbsSampler(dna, k, t, N):

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

        if CalculateScore(motifs) < CalculateScore(best_motifs):
            best_motifs = motifs
            best_score = CalculateScore(motifs)

        # Debug
        print "Round", str(j) + ", replacing motif", i, "with", new_motif, "to give",
        print "score:", CalculateScore(motifs), "Data:", motifs

    print "Best Score:", best_score

    return best_motifs

# End of GibbsSampler()

