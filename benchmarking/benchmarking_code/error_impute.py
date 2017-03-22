from __future__ import division
import numpy as np
from numpy import linalg as LA
from scipy.linalg import sqrtm, inv
import numpy as np
from sklearn.decomposition import PCA

class error(object):
    


    def RMSE(org,imputed,mask):
        return (((imputed[mask] - org[mask]) ** 2).mean())**(.5)
    
    def MSE(org,imputed,mask):
        return ((imputed[mask] - org[mask]) ** 2).mean()
    
    def forb(org,imputed,mask):
        return np.linalg.norm(imputed[mask] - org[mask],ord=2)/np.linalg.norm(org[mask],ord=2)
    
    def get_density(data):
        nonzeroscount=np.count_nonzero(data)
        sizel = data.shape
        totalentr=sizel[0]*sizel[1]
        return (nonzeroscount/totalentr*100)

    def get_pca_var(data):


    
        observed_table_sni, mapping = match(observed_table_sni, mappingdf)
        pca_model=PCA(n_components=3)
        X2=observed_table_sni.as_matrix()
        X_reduced2 = pca_model.fit_transform(X2)
        var_exp=pca_model.explained_variance_ratio_
        cum_var_exp = np.cumsum(pca_model.explained_variance_ratio_)

        return var_exp,cum_var_exp
