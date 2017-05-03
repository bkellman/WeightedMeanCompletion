# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:35:05 2017

@author: Benjamin Kellman
"""

from scipy.spatial.distance import pdist,squareform
import numpy as np

#' @export 
def sim_linear_func(D,c=1):
    1-(c*squareform(D))

def sim_exp_func(D,e=2,c=1):
    e**(1-(c*squareform(D)))

# d_i and d_j are distance matrixes
def interp_weightedMean(X,t=2,alpha=.5,D_i=None,D_j=None): #,sim_func_i=None,sim_func_j=None):
    if(alpha>1 or alpha<0):
        raise Exception('alpha must be between 0 and 1')
    if((D_i!=None and type(D_i)!=np.ndarray) or (D_j!=None and type(D_j)!=np.ndarray)):
        raise Exception('If specified, D_i and D_j must be a square distance matrix, class: np.ndarray')
    if(t<=0):
        raise Exception('t must be >=1')
    if(type(X)!=np.ndarray):
        raise Exception('X must be an np.ndarray')
    
    if(D_i==None):
        D_i = pdist(X,method='euclidean')
    if(D_j==None):
        D_j = pdist(X.transpose,method='euclidean')
    
    sim_func_i = sim_lin_func #(x,e=10) 
    sim_func_j = sim_exp_func #(x,e=100) 
    c_i = 10
    e_j = 100
  
    # normailize dist object
    D_i = D_i/D_i.max
    D_j = D_j/D_j.max
    # get similarity from distance objects
    W_i = sim_func_i(D_i,c=c_i)
    W_j = sim_func_j(D_j,e=e_j)
    # get normalization factors
    Omega_i = np.zeros(W_i.shape)
    diag_i = [1/sum(k) for k in range(W_i.shape[0])]
    Omega_i = np.fill_diagonal(Omega_i,diag_i)
    W_i_star = np.dot(Omega_i , W_i )

    Omega_j = np.zeros(W_j.shape)
    diag_j = [1/sum(k) for k in range(W_i.shape[1])]
    Omega_j = np.fill_diagonal(Omega_j,diag_j)
    W_j_star = np.dot(Omega_j , W_j )
  
    # smooth
    #X[is.na(X)]=0 if we end up with na's we will have to get rid of them
    while(t>0):
        t-=1
        X =  alpha * np.dot(W_i_star , X ) +  (1-alpha) * np.dot(X , W_j_star)
  return X
}
