# WeightedMeanCompletion

## Curated Phenotype and Phylogenic Database

- Both databases match in species represented with a total of 742 bacteria representing the type strain in the species.
- The Phenotype database has a total of 156 metadata factors flattened from the [DSMZ](https://bacdive.dsmz.de) api. 
- The Phylogenic database is based on cophenetic distance between the 742 species.  

### Collect Phenotype data

- Curated phenotype information was collected from bacdive [DSMZ](https://bacdive.dsmz.de) for type strain species. 

### Collect Phylogeny data

- A cophenetic distance matrix was curated based on the available species from [Phylophlan](https://huttenhower.sph.harvard.edu/phylophlan) 
- The benefit of matching species to Phylophlan is two fold. First, the Phylophlan included taxonomy represent the best classified species in each genus. While mainting a good phylogenic spread of all bacterial groups. Second, this method allows for "new" bacterial species to be accurately inserted into the distance matrix based on the Phylophlan genome insertion method. 

### Coding Tasks

- perform standard matrix-completion/interpolation: EM, regulariztion, nuclear norm...
- machine learning approach to matrix-completion/interpolation
- Weighted Mean Completion approach to matrix-completion/interpolation

# Paper
https://www.overleaf.com/6466536zjchwf#/21826208/
