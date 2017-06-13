# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:35:05 2017
@author: Benjamin Kellman
"""

from scipy.spatial.distance import pdist,squareform
import numpy as np
from sklearn import preprocessing
from itertools import combinations

#######
#from rpy2.robjects.packages import importr
#import rpy2.robjects as ro
#import pandas.rpy.common as com
#import rpy2.robjects.numpy2ri
#from rpy2.robjects import pandas2ri
#import rpy2.rlike.container as rpc
#pandas2ri.activate()
#rpy2.robjects.numpy2ri.activate()
#load package
#WMI = importr('WeightedMeanInterpolation')
from rpy2.robjects.packages import importr
stats = importr('stats')
base = importr('base')
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

########


#' @export 
def sim_linear_func(D,c=1):
    return 1-(c*D)

def sim_exp_func(D,e=2,c=1):
    return e**(1-(c*D))

# d_i and d_j are distance matrixes
def interp_weightedMean(X,t=2,alpha=.5,D_i=None,D_j=None): #,sim_func_i=None,sim_func_j=None):
    
    if alpha>1 or alpha<0:
        raise Exception('alpha must be between 0 and 1')
    if D_i!=None and type(D_i)!=np.array or D_j!=None and type(D_j)!=np.array:
        raise Exception('If specified, D_i and D_j must be a square distance matrix, class: np.ndarray')
    if t<=0:
        raise Exception('t must be >=1')
    #if type(X)!=np.array:
    #      raise Exception('X must be an np.ndarray')
    print('X')
    print(X[1])
    if D_i==None:
        #D_i = pdist(X,metric='euclidean')
        D_i = np.asarray(stats.dist(X,method='euclidean',upper=True,diag=True))
    if D_j==None:
        #D_j = pdist(X.T,metric='euclidean')
        D_j = np.asarray(stats.dist(X.T,method='euclidean',upper=True,diag=True))
    #print('dist')
    #print(D_i)
    #print(type(D_i))
    #print(np.asarray(D_i))
    #sim_func_i = sim_lin_func() #(x,e=10)
    #sim_func_j = sim_exp_func() #(x,e=100)
    c_i = 10
    e_j = 100
    # normailize dist object
    min_max_scaler = preprocessing.MinMaxScaler()
    D_i = min_max_scaler.fit_transform(D_i)
    D_j = min_max_scaler.fit_transform(D_j)
    #D_i = D_i/np.max(D_i)
    #D_j = D_j/np.max(D_j)
    # get similarity from distance objects
    W_i = sim_linear_func(D_i,c=c_i)
    W_j = sim_exp_func(D_j,e=e_j)
    # get normalization factors
    Omega_i = np.zeros(W_i.shape)
    diag_i = [1/W_i[k].sum() for k in range(0,W_i.shape[0])]
    Omega_i = np.fill_diagonal(Omega_i,diag_i)
    print(W_i[1])
    W_i_star = np.dot(Omega_i , W_i )

    Omega_j = np.zeros(W_j.shape)
    diag_j = [1/W_j[:,k].sum() for k in range(W_j.shape[1])]
    Omega_j = np.fill_diagonal(Omega_j,diag_j)
    W_j_star = np.dot(Omega_j , W_j )
  
    # smooth
    #X[is.na(X)]=0 if we end up with na's we will have to get rid of them
    while(t>0):
        print("Hello")
        t-=1
        X =  alpha * np.dot(W_i_star , X ) +  (1-alpha) * np.dot(X , W_j_star)
    return X
