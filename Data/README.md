## Data Files

- [SNP_positions.csv](SNP_positions.csv): Chromosome and position information of each TagPair as aligned to the Red deer reference genome (Bana et al. 2018)
- [Test.csv](Test.csv): Phenotypic data (Sex and Bread) for the Test dataset used in Bilton et al. (2018)
- [Train.csv](Train.csv): Phenotypic data (Sex and Bread) for the Training dataset used in Bilton et al. (2018)
- [deer_count.csv](deer_count.csv): Allele counts for the SNP-linked SNPs for the Training and Test datasets. The data is in TagDigger format such that every even column represents the counts for the reference allele of the named SNP and every odd column (except the first column which gives the sample name) represents the counts for the alternate allele of the named SNP.
- [deer_tagdigger.csv](deer_tagdigger.csv): Reference file matching the SNP names in the *deer_count.csv* with the TagPair names in the *SNP_positions.csv* file.

### References

Bana N.&#193;., Nyiri A., Nagy J., Frank K., Nagy T., St&#233;ger V., Schiller M., Lakatos P., Sug&#225;r L., Horn P., Barta E. & Orosz L. (2018) The red deer Cervus elaphus genome CerEla1.0: sequencing, annotating, genes, and chromosomes. *Molecular Genetics and Genomics* *293*, 665-84.

Bilton, T.P., Chappell, A.J., Clarke, S.M., Brauning, B., Dodds, K.G., McEwan, J.C. \& Rowe, S.J. (2018). Using genotyping-by-sequencing to predict gender in animals. Unpublished Manuscript

