# WeightedMeanCompletion

## Quick Start

```{R}
# installation
library(devtools)
install_github('bkellman/WeightedMeanCompletion')
library(WeightedMeanCompletion)

# interpolation in mtcars
X = data.matrix(mtcars)
X[sample(prod(dim(X)),20)] = NA
# default correlation based interpolation for rows (i) and columns (j)
X_out1 = interp_weightedMean(X,t=5)
# default correlation-based interpolation for columns (j), integer distance informed interpolation with exponential spreading function in rows (i)
X_out2 = interp_weightedMean(X,t=10,D_i=dist_integer(X),sim_func_i=sim_exp_func)
# integer distance informed interpolation with exponential spreading function in rows (i) and columns (j)
X_ou32 = interp_weightedMean(X,t=10,
  D_i=dist_integer(X),sim_func_i=sim_exp_func,
  D_j=dist_integer(X),sim_func_j=sim_exp_func)
```
For more information on parameter optimization, see the notes at the end of each function in [R/benchmarking.R](https://github.com/bkellman/WeightedMeanCompletion/blob/master/R/benchmarking.R)

## Curated Phenotype and Phylogenic Database Summary 

- Both databases match in species represented with a total of 742 bacteria representing the type strain in the species (50% match between databases).
- The Phenotype database has a total of 156 metadata factors flattened from the [Bacdive DSMZ api](https://bacdive.dsmz.de). 
- The Phylogenic database is based on cophenetic distance between the 742 species.  

### Phenotype data

- Curated phenotype information was collected from the [Bacdive DSMZ api](https://bacdive.dsmz.de) for type strain species. 

### Phylogeny data

- A cophenetic distance matrix was curated based on the available species from [Phylophlan](https://huttenhower.sph.harvard.edu/phylophlan) 
- The benefit of matching species to Phylophlan is two fold. First, the Phylophlan included taxonomy represent the best classified species in each genus. While mainting a good phylogenic spread of all bacterial groups. Second, this method allows for "new" bacterial species to be accurately inserted into the distance matrix based on the Phylophlan genome insertion method. 


## Benchmarking

### Benchmarking code
- base_impute: simple imputation methods such as mean,median,common fill, and zeros
- error_impute: can calculate RMSE,forbenious error and a function for measuring the explained varance of dataset through PCA
- benchmarking_methods: method that takes a dataset and runs 8 different methods for imputation for benchmarking 
** TODO: (must add our method to the methods benchmarked) 

### Benchmarking Dataset
- TODO: find Common imputation benchmarking dataset 
- Benchmarking performed on our dataset


## TODO

### Coding Tasks

- perform standard matrix-completion/interpolation: EM, regulariztion, nuclear norm...
- machine learning approach to matrix-completion/interpolation
- Weighted Mean Completion approach to matrix-completion/interpolation

# Paper
https://www.overleaf.com/6466536zjchwf#/21826208/
