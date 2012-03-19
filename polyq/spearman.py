#!/usr/bin/env python3

import numpy, os, scipy
from scipy.stats import spearmanr

os.chdir('C:/Users/user/Documents/holt polyq')

#Load csv
polyq = numpy.genfromtxt(
    'polyq_matched_x.csv', delimiter=',', skip_header=1, usecols=range(1,5))
warringer = numpy.genfromtxt('warringer_matched_x.csv',
                             delimiter=',', skip_header=1, usecols=range(1,5))


# Take only those counts where sd != 0
interestingQ = polyq[scipy.std(polyq, 1) != 0,]


result = numpy.zeros((len(interestingQ), len(warringer)))

# each gene's correlations for "fitness i" with qrow (polyq count)
def sp_corr_wrapper(qrow, i):
    result = spearmanr(qrow, warringer[i,])
    return result[0]

for i in range(len(warringer)):
    result[:, i] = numpy.apply_along_axis(sp_corr_wrapper, 1, interestingQ, i)

numpy.savetxt("correlations.csv", result, fmt='%6.4f', delimiter=',')

