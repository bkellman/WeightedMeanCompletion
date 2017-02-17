#low rank methods (new)
from wpca import WPCA, EMPCA
from fancyimpute import BiScaler, KNN, NuclearNormMinimization, SoftImpute, IterativeSVD, MICE, MatrixFactorization
from sklearn.decomposition import PCA
#other
from __future__ import division
import numpy as np
import pandas as pd


class benchmarking_methods(object):

    '''''
    This contains code for the benchmarking of imputation methods
    
    '''''

    def benchmark_complete(data,ending_density=.02,step):

    '''
        Input: Data array to benchmark on, the ending density to return results, the step bteween density imputation
        
        Output: Dataframe of output density and RMSE for each method with respect to each input density
        
        
        '''
        ############################  makes most sense (min only) #######################
        
        # removes min value that is greater than zero (checks density) in each iteration randomly chosen
        
        #density range to run
        nonzeroscount=np.count_nonzero(data)
        sizel = data.shape
        totalentr=sizel[0]*sizel[1]
        end=0.02 # final density to test
        begin=(nonzeroscount/totalentr) # Begning density of matrix given
        step=.01 # step of density
        
        
        #intialize lists to store
        density_in=[]
        RMSE_empca_scores=[]
        RMSE_wpca_scores=[]
        RMSE_sfi_scores=[]
        RMSE_siv_scores=[]
        RMSE_sni_scores=[]
        RMSE_smi_scores=[]
        RMSE_szi_scores=[]
        Density_empca=[]
        Density_wpca=[]
        Density_sfi=[]
        Density_siv=[]
        Density_sni=[]
        Density_smi=[]
        Density_szi=[]
        
        #radnomly remove values from known matrix and try to impute them
        
        for d in reversed(np.arange(end,begin,step)):
            otum=data.T.copy()
            
            #begin density check
            nonzeroscount=np.count_nonzero(otum)
            sizel = otum.shape
            totalentr=sizel[0]*sizel[1]
            
            while np.float64((nonzeroscount/totalentr)) > d:
                #remove a min frequency OTU and then check density
                j=randint(0,len(otum[:][:])-1)
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
            
            print("Running EMPCA")
            EMPCAi = EMPCA(n_components=3).fit_reconstruct(otum.copy(),weight)
            print("Running WPCA")
            WPCAi = WPCA(n_components=3).fit_reconstruct(otum.copy(),weight)

            # fancy impute and zeros
            print("Nuclear Norm")
            sni=NuclearNormMinimization(min_value=(np.amin(otum2)),max_value=(np.amax(otum2))).complete(otum2.copy())
            print("Running Soft Impute")
            sfi=SoftImpute(shrinkage_value=None,convergence_threshold=0.00001,max_iters=1000,max_rank=min(otum2.shape),n_power_iterations=1,init_fill_method="zero",min_value=(np.amin(otum2)),max_value=(np.amax(otum2)),normalizer=None,verbose=False).complete(otum2.copy())
            print("Running Iterative SVD")
            siv=IterativeSVD(rank=(min(otum2.shape)-1),convergence_threshold=0.00001,max_iters=1000,gradual_rank_increase=True,svd_algorithm="arpack",init_fill_method="zero",min_value=(np.amin(otum2)),max_value=(np.amax(otum2)),verbose=False).complete(otum2.copy())
            print("Running Matrix Factorization")
            smi=MatrixFactorization(rank=(min(otum2.shape)-1),initializer=np.random.randn,learning_rate=0.01,patience=5,l1_penalty=0.05,l2_penalty=0.05,min_improvement=0.01,max_gradient_norm=5,optimization_algorithm="adam",min_value=(np.amin(otum2)),max_value=(np.amax(otum2)),verbose=False).complete(otum2.copy())
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
            Density_sni.append(error.get_density(sni))
            Density_smi.append(error.get_density(smi))
            Density_szi.append(error.get_density(szi))
            
            # RMSE of imputed values
            missing_mask = np.isnan(otum2.T) # masking to only check RMSE between values imputed and values removed
            RMSE_empca_scores.append(error.RMSE(data,EMPCAi,missing_mask))
            RMSE_wpca_scores.append(error.RMSE(data,WPCAi,missing_mask))
            RMSE_sfi_scores.append(error.RMSE(data,sfi.T,missing_mask))
            RMSE_siv_scores.append(error.RMSE(data,siv.T,missing_mask))
            RMSE_sni_scores.append(error.RMSE(data,sni.T,missing_mask))
            RMSE_smi_scores.append(error.RMSE(data,smi.T,missing_mask))
            RMSE_szi_scores.append(error.RMSE(data,szi.T,missing_mask))
        
        
        RMSEmapping = pd.DataFrame({'Density': list(map(int, density_in)),'EMPCA': RMSE_empca_scores,'Matrix Factorization': RMSE_smi_scores,'WPCA': RMSE_wpca_scores,'Soft Impute': RMSE_sfi_scores,'Iterative SVD': RMSE_siv_scores,'Nuclear Norm Minimization': RMSE_sni_scores,'Zeros Replace Unknown': RMSE_szi_scores})
        RMSEmapping.set_index(['Density'], inplace=True)
        Out_density = pd.DataFrame({'density': list(map(int, density_in)),'EMPCA': Density_empca,'Matrix Factorization': Density_smi,'WPCA': Density_wpca,'Soft Impute': Density_sfi,'Iterative SVD': Density_siv,'Nuclear Norm Minimization': Density_sni,'Zeros Replace Unknown': Density_szi})
        Out_density.set_index(['density'], inplace=True)
        
        return Out_density, RMSEmapping

