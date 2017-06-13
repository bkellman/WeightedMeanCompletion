from __future__ import division

#low rank methods (new)
from wpca import WPCA, EMPCA
from fancyimpute import BiScaler, KNN, NuclearNormMinimization, SoftImpute, IterativeSVD, MICE, MatrixFactorization
from sklearn.decomposition import PCA
from sklearn import preprocessing
#other
import numpy as np
import pandas as pd
#local
from base_impute import base
from error_impute import error
import core

class benchmarking_methods(object):

    '''''
    This contains code for the benchmarking of imputation methods
    
    '''''

    def benchmark_complete(data,ending_density=.02,step=.01):

        '''
        Input: Data array to benchmark on, the ending density to return results, the step bteween density imputation
        
        Output: Dataframe of output density and RMSE for each method with respect to each input density
        
        
        '''
        # removes min value that is greater than zero (checks density) in each iteration randomly chosen
        #density range to run
        nonzeroscount=np.count_nonzero(data)
        sizel = data.shape
        totalentr=sizel[0]*sizel[1]
        end=0.02 # final density to test
        begin=(nonzeroscount/totalentr) # Begning density of matrix given
        #step=.01 # step of density
        
        
        #intialize lists to store
        density_in=[]
        RMSE_empca_scores=[]
        RMSE_wpca_scores=[]
        RMSE_sfi_scores=[]
        RMSE_siv_scores=[]
        RMSE_szi_scores=[]
        RMSE_wmiC_scores=[]
        RMSE_wmiP_scores=[]
        Density_empca=[]
        Density_wpca=[]
        Density_sfi=[]
        Density_siv=[]
        Density_szi=[]
        Density_wmiC=[]
        Density_wmiP=[]
        
        #radnomly remove values from known matrix and try to impute them
        
        for d in reversed(np.arange(end,begin,step)):
            otum=data.T.copy()
            
            #begin density check
            nonzeroscount=np.count_nonzero(otum)
            sizel = otum.shape
            totalentr=sizel[0]*sizel[1]
            
            while np.float64((nonzeroscount/totalentr)) > d:
                #remove a min frequency OTU and then check density
                j=np.random.randint(0,len(otum[:][:])-1)
                #make sure row is not all zero (all zero row causes singular matrix)
                if sum(list(otum[j][:])) < 1:
                    continue
                m = min(i for i in list(otum[j][:]) if i > 0)
                #make sure removing value will not result in zero row
                if sum(list(otum[j][:])) == m:
                    continue
                otum[j][list(otum[j][:]).index(m)]=0
                #check denstiy to break
                nonzeroscount=float(np.count_nonzero(otum))
                sizel = otum.shape
                totalentr=float(sizel[0])*float(sizel[1])
            
            
            # coherce float of the unknown and print new density
            print("Data table of %f generated"%d)
            otum=otum.T.astype(np.float64)
            
            # make zero unknown for fancy impute, avoid singular matrix by taking transpose
            otum2=otum.T.copy()
            otum2=otum2.astype(np.float64)
            otum2[otum2 == 0] = np.nan #make unknown nan
            
            #WPCA and EMPCA
            
            #build wieghted matrix
            weight = otum.copy()
            for i in range(len(otum2.T)):
                for j in range(len(otum2.T[i])):
                    if otum2.T[i][j]==0:
                        weight[i][j]=1
                    else:
                        weight[i][j]=1000
        
        
            print("Weighted Mean Interpolation without phylo-distance")
            otum2=np.asarray(otum2)
            wmiC=core.interp_weightedMean(otum2.copy())
            print("Weighted Mean Interpolation with phylo-distance")
            phylo = pd.read_csv('data/Matched_Pheno_and_Phylo_Data/matched_phylo.csv/matched_phylo.csv')
            wmiP=core.interp_weightedMean(otum2.copy(),D_i=phylo.as_matrix())
        
            print("Running EMPCA")
            EMPCAi = EMPCA(n_components=3).fit_reconstruct(otum.copy(),weight)
            print("Running WPCA")
            WPCAi = WPCA(n_components=3).fit_reconstruct(otum.copy(),weight)

            # fancy impute and zeros
            print("Running Soft Impute")
            sfi=SoftImpute(verbose=False).complete(otum2.copy())
            print("Running Iterative SVD")
            siv=IterativeSVD(verbose=False).complete(otum2.copy())
            print("Imputing by filling with zeros for base comparison")
            szi=base.zeros(otum2.copy())

            # save the results
    
            #density in (after removed values)
            density_in.append(error.get_density(otum))
            
            # density imputed
            Density_empca.append(error.get_density(EMPCAi))
            Density_wpca.append(error.get_density(WPCAi))
            Density_sfi.append(error.get_density(sfi))
            Density_siv.append(error.get_density(siv))
            Density_szi.append(error.get_density(szi))
            Density_wmiC.append(error.get_density(wmiC))
            Density_wmiP.append(error.get_density(wmiP))
            
            # RMSE of imputed values
            missing_mask = np.isnan(otum2.T) # masking to only check RMSE between values imputed and values removed
            RMSE_empca_scores.append(error.RMSE(data,EMPCAi,missing_mask))
            RMSE_wpca_scores.append(error.RMSE(data,WPCAi,missing_mask))
            RMSE_sfi_scores.append(error.RMSE(data,sfi.T,missing_mask))
            RMSE_siv_scores.append(error.RMSE(data,siv.T,missing_mask))
            RMSE_szi_scores.append(error.RMSE(data,szi.T,missing_mask))
            RMSE_wmiC_scores.append(error.RMSE(data,wmiC.T,missing_mask))
            RMSE_wmiP_scores.append(error.RMSE(data,wmiP.T,missing_mask))
        
        
        RMSEmapping = pd.DataFrame({'Density': list(map(int, density_in)),'EMPCA': RMSE_empca_scores,'WPCA': RMSE_wpca_scores,'Soft Impute': RMSE_sfi_scores,'Iterative SVD': RMSE_siv_scores,'Zeros Replace Unknown': RMSE_szi_scores,'Weighted-Mean Interpolation Correlation':RMSE_wmiC_scores,'Weighted-Mean Interpolation Phylo':RMSE_wmiP_scores})
        RMSEmapping.set_index(['Density'], inplace=True)
        Out_density = pd.DataFrame({'density': list(map(int, density_in)),'EMPCA': Density_empca,'WPCA': Density_wpca,'Soft Impute': Density_sfi,'Iterative SVD': Density_siv,'Zeros Replace Unknown': Density_szi,'Weighted-Mean Interpolation Correlation':Density_wmiC,'Weighted-Mean Interpolation Phylo':Density_wmiP})
        Out_density.set_index(['density'], inplace=True)
        
        return Out_density, RMSEmapping

