#!/usr/bin/env python

import pandas as pd
from benchmarking_methods_no_nuclear import benchmarking_methods
from sklearn import preprocessing

pheno_data=pd.read_table('/Users/cameronmartino/bin/WeightedMeanCompletion/data/Matched_Pheno_and_Phylo_Data/matched_pheno.csv/matched_pheno.csv',index_col=0,low_memory=False)
pheno_data[list(pheno_data.columns.values)] = pheno_data[list(pheno_data.columns.values)].astype(str)
pheno_data=pheno_data.dropna(axis=1,how='all')
pheno_data=pheno_data.dropna(axis=0,how='all')

pheno_data=pheno_data.apply(preprocessing.LabelEncoder().fit_transform)
Out_density, RMSEmapping = benchmarking_methods.benchmark_complete(pheno_data.as_matrix())

Out_density.to_csv('/Users/cameronmartino/bin/WeightedMeanCompletion/benchmarking/benchmarking_code/Out_density.csv', sep='\t')
RMSEmapping.to_csv('/Users/cameronmartino/bin/WeightedMeanCompletion/benchmarking/benchmarking_code/RMSEmapping.csv', sep='\t')
