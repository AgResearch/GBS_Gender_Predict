## R Scripts

R code for analyzing genotyping-by-sequencing data and predict the gender of red deer. There are two parts to this program. 

1. [find_SNPs_function.R](find_SNPs_function) loads in the function `finding_SNPs()`. The purpose of this function is to find the tagPairs that have aligned to the sex chromosomes of the reference genome (for the species). It then finds the matching marker from the SNP catalogue and checks the behaviour of the SNPs in KGD with a set of sex verified samples. It removes SNPs that do not meet certain criteria. This code only needs to be run once each time a new SNP catalogue is made. The argument for this function are:
    - `SNP_positions`: Filename to a csv file with the chromosome and position information for each TagPair.
    - `deer_tagdigger`: Filename of a csv filegiving the matching TagPairs names from the `SNP_positions` file with the SNP names used in the `deer_count` file.
    - `deer_count`: Filename for the allele count file in TagDigger format.
    - `test_samples`: Filename of a csv containing the phenotypic data for the training datasets to verify the SNPs.
    - `output_file`: Path to a directory in which the results are to be written to.

2. [predict_gender](predict_gender.R) loads in the function `predict_gender()` which performs the gender prediction. The arguments of this function are:
    - `sampleID`: Filename of a csv file containing the ID's and phenotypic data for the dataset to predict gender for.
    - `deer_results`: Filename of a csv file containing the read counts for the dataset subsetted to contain only the sex-linked markers only (produced from the `finding_SNPs()` function).
    - `sex_markers`: Filename of a csv file containing the information of the SNPs determined to be on the Sex chromosome (produced from the `finding_SNPs()` function)
    - `filtered_markers`: Filename of a csv file containing the column names of the of markers to remove (produced from the `finding_SNPs()` function).
    - `output_file`: Path to a directory in which the results are to be written to.
    - `withPlotly`: Logical argument specifying whether to produce interactive plots using plotly (if `TRUE`).
    - `group`: A character string giving the column name of the `sampleID` file to use as a grouping variable.

The third script, [GBS-Chip-Gmatrix.R](GBS-Chip-Gmatrix.R), is a copy of the [KGD](https://github.com/AgResearch/KGD) program (v0.2) and is reproduced (with a few minor changes) for ease of running the code.