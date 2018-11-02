
## Filter training data set to remove erroreous and ficitictious SNPs
## and run through KGD
source("../Scripts/find_SNPs_function.R")  #code to source the predict gender function 

SNP_positions <- "../Data/SNP_positions.csv"
deer_tagdigger <- "../Data/deer_tagdigger.csv"
deer_count <- "../Data/deer_count.csv"
sampleID <- "../Data/Train.csv"

finding_SNPs(SNP_positions,deer_tagdigger,deer_count,sampleID,output="./")


##### Predict the gender for training and test data sets
source("../Scripts/Gender_Prediction_Function.R")  #code to source the predict gender function 

deer_results <- "deerAll20170301_counts_filtered_sexmarkers_only.csv"
sex_markers <- "sex_markers.csv"
filtered_markers <- "markers_to_remove.csv"

## Training data set
sampleID <- "../Data/Train.csv"
predict_gender(sampleID,deer_results,sex_markers,filtered_markers,output_file="Train", withPlotly=T)

## Red deer test data set
sampleID <- "../Data/Test.csv"
predict_gender(sampleID,deer_results,sex_markers,filtered_markers,output_file="Test", withPlotly=T, group="Breed")

