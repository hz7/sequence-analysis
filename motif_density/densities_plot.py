'''
Created on Feb 24, 2012

@author: user
'''
import h5py, numpy
import matplotlib.pyplot as plt

top = [('YMR273C', 0.28000000000000003), ('YIL154C', 0.20000000000000001), ('YMR070W', 0.20000000000000001), ('YEL043W', 0.20000000000000001), ('YKR075C', 0.20000000000000001), ('YOR247W', 0.20000000000000001), ('YDL055C', 0.16), ('YJR117W', 0.16), ('YGL173C', 0.16), ('YOL155C', 0.16), ('YAR042W', 0.16), ('YIR036C', 0.16), ('YOR070C', 0.16), ('YPR181C', 0.16), ('YOL161C', 0.16), ('YCR079W', 0.16), ('YGR227W', 0.16), ('YMR296C', 0.16), ('YLL066W-B', 0.16), ('YML059C', 0.16)]

f = h5py.File('densities.h5', 'r')

for pair in top:
    key = pair[0]
    data0 = f[key]['0'][0]
    data1 = f[key]['1'][0]
    data2 = f[key]['2'][0]
    plt.figure()
    plt.plot(numpy.arange(len(data0)), data0,
             numpy.arange(len(data1)), data1,
             numpy.arange(len(data2)), data2)
    plt.title(key)
    plt.legend(('shift 0', 'shift 1', 'shift 2'))
    plt.savefig('density_' + key + '.png', format='png')
    