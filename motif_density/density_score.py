'''
pick some genes where densities is high (scored by max density)
'''

import h5py, numpy

f = h5py.File('densities.h5', 'r')

score_table = {}

for key in f.keys():
    data1 = f[key]['1']
    data2 = f[key]['2']
    if len(data1) > 0:
        score1 = numpy.max(data1[...])
    else:
        score1 = -1
    if len(data2) > 0:
        score2 = numpy.max(data2[...])
    else:
        score2 = -1
    score_table[key] = max(score1, score2)

sorted_table = sorted(score_table.items(), key=lambda pair: pair[1], reverse=True)
print(sorted_table[:40])
