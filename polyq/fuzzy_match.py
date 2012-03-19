#!/usr/bin/env python3

"""
Try to match common species.
Might have to use a database.
For two ids, idA and idB:
idA - {Asynonym1, Asynonym2, ...}
idB - {Bsynonym1, ...}
try to match synonyms?

seems useful:
http://www.uniprot.org/taxonomy/764097
"""

import csv, difflib, os, re

os.chdir('C:/Users/user/Documents/holt polyq')

with open('warringer.csv', 'r') as f:
	warringerRead = csv.reader(f)
	warringerHeader = next(warringerRead)

with open('polyq.csv', 'r') as f:
	polyqRead = csv.reader(f)
	polyqHeader = next(polyqRead)

warringerStrains = warringerHeader[4:]
polyqStrains = [re.sub('_.+$', '', x) for x in polyqHeader[1:]]

for Qstrain in polyqStrains:
    matches = difflib.get_close_matches(Qstrain, warringerStrains, cutoff=0.5)
    if len(matches) > 0:
        print([Qstrain, matches])

