'''
pick some genes where min. distance is interesting -- whatever that means
'''

import h5py, numpy

def score(distances):
    if len(distances) == 0:
        return float("inf")
    else:
        nDist = numpy.size(distances)
        min5 = numpy.sort(distances)[:min(5, nDist)]
        return numpy.average(min5) + max(5 - nDist, 0)

f = h5py.File('distance_results.h5', 'r')

score_table = {}

for key in f.keys():
    data1 = f[key]['1']
    data2 = f[key]['2']
    if len(data1) > 0:
        score1 = score(data1[...])
    else:
        score1 = float("inf")
    if len(data2) > 0:
        score2 = score(data2[...])
    else:
        score2 = float("inf")
    score_table[key] = min(score1, score2)

sorted_table = sorted(score_table.items(), key=lambda pair: pair[1])
print(sorted_table[:40])
