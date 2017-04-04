#!/usr/bin/env python

import sys
import os
import time

import utils

def ReadGenome(genome_file_name):

    genome = ''
    with open(genome_file_name, 'r') as genome_fh:
        for line in genome_fh.readlines():
            if not line.startswith('>'):
                genome += line.strip()

    return genome

# End of ReadGenome()

def FindDnaABox(genome):

    minimum_skews = utils.FindMinimumSkews(genome)

    for minimum_skew in minimum_skews:
        print("Minimum Skew:", minimum_skew)

    # Minimum skews for Salmonells enterica are 3764856 and 3764858.  Ideally
    # we should check all minimi, or if some are close to each other, combine
    # them into one mininum, but since we're only doing S. enterica for now,
    # and since it's two minimi are so close to each other, I'm just going to
    # pick the first one and go from there.

    minimum_skew = minimum_skews[0]

    # Find most frequent 9-mers with 1 mismatch in area +-500 around
    # predicted oriC.
    region_length = 1000
    region_start = minimum_skew - int(region_length/2)
    region_end = region_start + region_length
    region = genome[region_start:region_end]
    patterns = utils.FindMostFrequentWithMismatchesAndReverseComplementWithLocci(9, 1, region)
    for pattern, locations in patterns.items():
        print(pattern, "found at locci:", [x + region_start for x in locations])

    return 'TBD'

# End of FindDnaABox()


if __name__ == "__main__":

    start_time = time.time()
    genome_file_name = sys.argv[1]

    genome = ReadGenome(genome_file_name)

    DnaABox = FindDnaABox(genome)
    end_time = time.time()

    print(DnaABox)
    print("Elapsed Time:", end_time - start_time)
