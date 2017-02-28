#!/usr/bin/env python

import sys
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


