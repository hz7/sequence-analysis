'''
plot densities of top min. distances
'''
import h5py, numpy
import matplotlib.pyplot as plt

top = [('YHR006W', 2.2000000000000002), ('YIL082W-A', 2.3999999999999999), ('YJR147W', 2.6000000000000001), ('YFL021W', 2.6000000000000001), ('YOR116C', 3.0), ('YPL054W', 3.0), ('YKL048C', 3.0), ('YGR109W-B', 3.2000000000000002), ('YIL169C', 3.3999999999999999), ('YJL141C', 3.6000000000000001), ('YGL066W', 3.666666666666667), ('YBR038W', 4.0), ('YEL047C', 4.0), ('YIL002C', 4.0), ('YDR489W', 4.0), ('YAR042W', 4.0), ('YPL122C', 4.0), ('YCR097W', 4.0), ('YGR024C', 4.0), ('YMR052W', 4.0), ('YHR030C', 4.0), ('YLR003C', 4.0), ('YLR154C-G', 4.0), ('YNR069C', 4.0), ('YGL100W', 4.0), ('YBR176W', 4.0), ('YML034W', 4.0), ('YPR005C', 4.0), ('YLR098C', 4.0), ('YIR019C', 4.0), ('YKL185W', 4.0), ('YDR501W', 4.0), ('YPL020C', 4.0), ('YMR084W', 4.0), ('YPL104W', 4.0), ('YLL066W-B', 4.0), ('YPL236C', 4.0), ('YGR225W', 4.0), ('YGR142W', 4.0), ('YLR201C', 4.0)]

f = h5py.File('densities.h5', 'r')

i = 1
for pair in top:
    key = pair[0]
    score = pair[1]
    data0 = f[key]['0'][0]
    data1 = f[key]['1'][0]
    data2 = f[key]['2'][0]
    plt.figure()
    plt.plot(numpy.arange(len(data0)), data0,
             numpy.arange(len(data1)), data1,
             numpy.arange(len(data2)), data2)
    plt.title('gene:{0}, rank:{1}, score:{2}'.format(key, i, score))
    plt.legend(('shift 0', 'shift 1', 'shift 2'))
    plt.savefig('distance_rank{0}_{1}.{2}'.format(i, key, 'png'), format='png')
    i += 1