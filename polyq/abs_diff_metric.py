'''
SUM of all absolute differences from one species to the next. e.g. If the table looked like this:

Density:         1, 4, 6, N/A, 6, 3, 10
Abs Dif:              3, 2, 0, 3, 7        
Sum of abs dif = 15

Where N/A is a missing value - here you just skip that position and go to the next.
'''

import numpy, pandas
from pandas import DataFrame

df = pandas.read_csv('c:/holt/polyq/data/out_pandas.csv', index_col=0)

result = DataFrame({'abs_diff' : [numpy.sum(numpy.abs(numpy.diff(row.dropna().values)))
                    for rowIndex, row in df.iterrows()]}, index=df.index)

result.to_csv('c:/holt/polyq/data/out_pandas_diffmetric.csv')
