#Created on 7/17/14
__author__ = 'Juan A. Ugalde'

from cogent import LoadSeqs
import itertools
from collections import defaultdict

aln = LoadSeqs("Acido251.faa")

alignment_length = len(list(aln.iterPositions()))

global_comp = defaultdict(int)

no_gaps = defaultdict(int)

for position in list(aln.iterPositions()):
    combinations = list(itertools.combinations(position, 2))

    for entry, pair in zip(combinations, range(len(combinations))):

        if entry[0] == entry[1]:
            global_comp[pair] += 1

        else:
            global_comp[pair] += 0

identity_list = []

for count in global_comp.values():
    identity_list.append(float(count)/alignment_length)

average_global_identity = sum(identity_list) / len(identity_list) * 100

print average_global_identity


