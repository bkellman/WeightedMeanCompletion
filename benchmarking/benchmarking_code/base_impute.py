from __future__ import division

try:
    import numpy as np
except ImportError:
    print('Unable to import numpy.')
try:
    from sklearn.preprocessing import Imputer
except ImportError:
    print('Unable to import Imputer from sklearn.preprocessing.')

from sklearn.preprocessing import Imputer
import numpy as np

#rpy2
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
import rpy2.rlike.container as rpc
pandas2ri.activate()
rpy2.robjects.numpy2ri.activate()
#load package
WMI = importr('WeightedMeanInterpolation')

class base(object):
    

    
    #Zero fill
    def zeros(A):
        return np.nan_to_num(A)
    #return A
    
    #mean fill benchmark
    def mean_fill(A):
        imp =Imputer(missing_values='NaN', strategy='mean', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    #median fill benchmark
    def median_fill(A):
        imp =Imputer(missing_values='NaN', strategy='median', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    #most_frequent fill benchmark
    def most_frequent_fill(A):
        imp =Imputer(missing_values='NaN', strategy='most_frequent', axis=0, verbose=0, copy=True)
        imp.fit(A)
        return imp.transform(A)
    
    def xrange(x):
        return iter(range(x))
    
    # weighted mean interpolation wrapper
    def wmi_wrapper(X,t=2,alpha=.6,**kwargs):
        return WMI.interp_weightedMean(X,t,alpha,rpc.OrdDict(kwargs))

