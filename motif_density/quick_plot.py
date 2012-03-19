'''
Created on Mar 8, 2012

@author: user
'''
import h5py, numpy
import matplotlib.pyplot as plt

f = h5py.File('data/densities_50.h5', 'r')

gene = 'YBL039W-B' #gene to plot

data0 = f[gene]['0'][0]
data1 = f[gene]['1'][0]
data2 = f[gene]['2'][0]
plt.figure()
plt.plot(numpy.arange(len(data0)), data0,
         numpy.arange(len(data1)), data1,
         numpy.arange(len(data2)), data2)
plt.title('gene: {}'.format(gene))
plt.legend(('shift 0', 'shift 1', 'shift 2'))
plt.show()
