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


if __name__ == "__main__":

    genome_file_name = sys.argv[1]

    genome = ReadGenome(genome_file_name)

    skew_diagram = utils.Skew(genome)

    genome_name = os.path.splitext(os.path.basename(genome_file_name))[0]
    with open('results/skew_diagram_for_' + genome_name + '.txt', 'w') as r_file:
        r_file.write(','.join([str(x) for x in skew_diagram]))

