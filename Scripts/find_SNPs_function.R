
library(data.table)



#finding SNP function, needs four arguments 
#--------------------------------------------------------------

finding_SNPs<-function(SNP_positions,deer_tagdigger,deer_count,test_samples, output_file="./"){

#output_file<- "/dataset/Predicting_Gender_Deer/active/Code/files_input_output_for_prediction"

SNP_pos<- read.csv(file=SNP_positions, header=TRUE) #read in a SNP position file
chromosome<- SNP_pos[,"chrom"] #creates a vector of chromosomes
X_SNPs<-which(chromosome == "CM008041.1") #find index of X chromosomes, name off NCBI
Y_SNPs<-which(chromosome == "CM008042.1") #find index of Y chromosomes, name off NCBI
XY_SNPs<- c(X_SNPs, Y_SNPs) #combining vectors of the indicies of X and Y chromosomes 


#prints total number of SNPs found and total numbers found on X and Y chromosomes seperately 
#----------------------------------------------------------------

cat(length(XY_SNPs), "sex SNPs found\n")
cat(length(X_SNPs) ,"X SNPs found\n")
cat(length(Y_SNPs) ,"Y SNPs found\n")

#finds rows of the sex_SNPs
#----------------------------------------------------------------

sex_SNPs<-SNP_pos[XY_SNPs,] #subsetting SNP_positoins to just rows of sex chromosomes
write.csv(sex_SNPs, file= file.path(output_file,"sex_SNP_positions.csv"),row.names=FALSE) #producing the outfile of just the rows of the X and Y chromosomes 

#reads in the file that has the tag digger generated marker names and the coresponding TP number to match the marker name with the
#TP number generated when calling the SNPs
#----------------------------------------------------------------

tag_digger<-read.csv(file=deer_tagdigger, header=TRUE)

subset_tagdigger<-tag_digger[tag_digger$UNEAK %in% sex_SNPs$tagpair,] #matches UNEAK tagpair to tagpair of sex chromosomes in sex_SNP_pos file
sex_markers<-subset_tagdigger[,"Marker.name"] #makes vector of just the marker names
row_index<-(match(subset_tagdigger$UNEAK,sex_SNPs$tagpair)) #using these indicies found above to extract the appropiate chromosomes
chromosome_find<-sex_SNPs[row_index,"chrom"] #using these indicies found above to extract the appropiate chromosomes
position_find<-sex_SNPs[row_index,"pos"] #using the indicies to extract apprpoiate positions 
subset_tagdigger[, "chromosome"]<-chromosome_find #adding the list of chromosomes onto the subset_tagdigger data.frame 
subset_tagdigger[, "SNP_position"]<-position_find #adding the list of SNP positions onto the subset_tagdigger data.frame 

if(!length(XY_SNPs) == length(sex_markers)){
cat("Number of marker names found does not equal number of TP numbers found")
}

cat(length(sex_markers), "sex marker names found in tagdigger file\n") #should be equal to total number of sex SNPs in previous script 
write.csv(subset_tagdigger, file = file.path(output_file,"sex_markers.csv"),row.names=FALSE) #creates a file of the sex markers,chromosome they appear on and position on the chromosome 

#subsetting the large deer file down by columns of sex markers #only 
#-----------------------------------------------------------

vector_sex_markers<-subset_tagdigger[,"Marker.name"] #creates a vector of just marker names

deer_count_header<-fread(deer_count,nrows=0,data.table=FALSE) #file name of the counts/reads in only the header row i.e the name of the markers 
marker_name<-colnames(deer_count_header) #creates a vector of marker names 
split_marker_name<-strsplit(as.character(marker_name), "_") #splits the name by underscore
sex_marker_name<-lapply(split_marker_name, "[[",1) #takes the first part of the string split
vector_marker_name<-unlist(sex_marker_name, use.names=FALSE) #turns into a vector
vector_cols<-which(vector_marker_name %in% vector_sex_markers) #find sex markers index
deer_count_names<-c("V1",marker_name[vector_cols])#finds names of sex_markers, include first column 
cat("Sub-setting deer count file by column,takes about three minutes\n")
subset_deer<-fread(deer_count, select= deer_count_names,  showProgress=TRUE ) 
#cat("Writing out file\n")
#write.csv(subset_deer, file= file.path(output_file, "deerAll20170301_counts_sexmarkers_only.csv")) #writes out file with sex SNP in columns only


#subsetting the above by rows to only include the samples we want 
#----------------------------------------------------------------

samples<-read.csv(file=test_samples, header=TRUE) #file can be changed to read in any file of samples with verified gender 
vector_of_samples<-samples[,"SampleID"] #create a vector of just the sampleIDs
cat(length(vector_of_samples),"samples loaded\n") #print number of samples found


cat("Now sub-setting deer count by sample\n")
SampleID<-subset_deer$V1 # making vector of the SampleIDs of this file
split_sample_ID<-strsplit(as.character(SampleID), "_") #breaks sampleID up by splitting at the underscores 
Sample_number<- lapply(split_sample_ID, "[[",1) #takes the first part of the entire sampleID
vector_sampleID<-unlist(Sample_number, use.names=FALSE) #turns this above list into a vector 

vector_rows<-which(vector_sampleID  %in% vector_of_samples) #finds the sampleIDs in the large file by matching with the sampleIDs found from script 03
deer_count_sub<-subset_deer[vector_rows,1:ncol(subset_deer)] #only taking rows of the samples we want in count file 
write.csv(deer_count_sub,file=file.path(output_file,"deer_count_subsetted_completely.csv"),row.names=FALSE)
cat(length(vector_rows), "samples found in deer count file") 


#loading in code to run KGD
#----------------------------------------------------------------

gform <- "TagDigger"  ####Used for RA files

#gform <- "uneak"   ####Used for HapMap files

#genofile <- "/dataset/Predicting_Gender_Deer/active/Code/files_input_output_for_prediction/deer_count_subsetted_completely.csv"   ###Add location to Ra or HapMap file
genofile <- file.path(output_file,"deer_count_subsetted_completely.csv")
sampdepth.thresh <- 0.3

source("../Scripts/GBS-Chip-Gmatrix.R",local=TRUE)


#merging samples from the same individuals, code from Ken Dodds 
#----------------------------------------------------------------

#merge up genotypes from samples that are from the same animal------------------------------------------------------------------
ped0<-read.csv(file=test_samples, header=TRUE) #need a file that has information on animalID 
SampleID <- read.table(text=seqID,sep="_",fill=TRUE,stringsAsFactors=FALSE)[,1]
pedpos <- match(SampleID,ped0$SampleID)
indsubset <- which(!is.na(pedpos))
animalID <- ped0$AnimalID[pedpos[indsubset]] #the variable after $ needs to be changed to the column name that holds ID

system.time(mg <- mergeSamples(animalID) )
nind <- mg$nind
seqID <- mg$seqID
genon <- mg$genon
depth.orig <- mg$depth.orig
sampdepth <- mg$sampdepth
snpdepth <- mg$snpdepth
pg <- mg$pg
animalID <- mg$mergeIDs
depth <-depth.orig
depth[depth < 2] <-1.1 
SampleID <-read.table(text=seqID,sep="_",fill=TRUE,stringsAsFactors=FALSE)[,1]

#End Merge Animals
#------------------------------------------------------------

SNPs_to_remove<-character() #creates an empty vector that will hold the names of all the SNPs that are needed to be removed

#need to know positions of the markers and what chromosome the markers appear on
#--------------------------------------------------------------------------------------------
markers<-match(SNP_Names, subset_tagdigger$Marker.name) #creates a vector of row indices/of the markers we want out of the sex_marker file
sex_markers<-subset_tagdigger[markers,] #holds only the marker names on the X and Y chromosome

#finding index of X and Y markers
#---------------------------------------------------------------------------------------------
index_Y_SNPs<-which(sex_markers$chromosome == "CM008042.1") #holds indicies of columns in the genon/ matrix which have the Y SNPs


index_X_SNPs<-which(sex_markers$chromosome == "CM008041.1")#holds indicies of columns in genon matrix/ with the X SNPs 


#finding the verified gender from our test set
#----------------------------------------------------------------------------------------------
verified_samples<-read.csv(file=test_samples, header=TRUE)
animals<-match(animalID,verified_samples$AnimalID) #finds index of the animals we want from test set
samples_in_genon<-verified_samples[animals,] #only holds animalIDs of those included after the merge
gender<-samples_in_genon$Sex
female_rows<-which(gender=="F") #holds index of the rows in genon that have female samples
male_rows<-which(gender=="M") #holds index of rows in genon that have male samples 

#test of behaviour for Y SNPs
#females should not have recorded SNP readings for SNPs on the Y chromosome 
#a Y SNP where females have recorded SNPs do not behave how we would expect and should be removed
#---------------------------------------------------------------- 
proportion_non_missing<-apply(genon[female_rows,index_Y_SNPs],2, function(x) sum(!is.na(x))/length(female_rows)) #for each SNP working out how many reads a female has, should be zero
uninformative_Y_SNPs<-which(proportion_non_missing>0.05) #allows for errors in alignment and homologous regions, finds the index of markers to remove
marker_uninformative_Y<-colnames(genon[0,index_Y_SNPs[uninformative_Y_SNPs]]) #gets the names of the columns to remove 
SNPs_to_remove<-append(SNPs_to_remove, marker_uninformative_Y) #adds marker names to vector



#another test for behaviour of SNPs on the Y chromosome
#males should not have heterozygous SNPs on the Y chromosome, SNPs should not have a high proportion of 1's
#ideally Y SNPs recorded for males shoud have zero 1's in the genon matrix
#-----------------------------------------------------------------------------------------------------------
proportion_hetero_Y<-apply(genon[male_rows,index_Y_SNPs],2,function(x) sum(x=="1")/length(male_rows)) #for each SNP calculating proportion of heterozygous reads
hetero_Y_SNPs<-which(proportion_hetero_Y>0.5) #removing SNPs where a fifth of the male samples are heterozygous, this may need to be altered with subsequent testing
marker_hetero_Y<-colnames(genon[0,index_Y_SNPs[hetero_Y_SNPs]])#gets the names of the columns to add vector
SNPs_to_remove<-append(SNPs_to_remove,marker_hetero_Y) #adds marker name to vector




#test to remove X SNPs that occur in the PAR (pseudo_autosomal region)
#will look at SNP positions and remove SNPs that occur on the X chromosome beyond 170,000,000 (estimated PAR boundary)
#------------------------------------------------------------------------------------------------------------
SNP_positions<-sex_markers$SNP_position  #creates a vector of the SNP positions 
SNPs_in_PAR<-which(SNP_positions>=170000000) #finds index of SNPs that occur beyond the threshold 
markers_in_PAR<-colnames(genon[0,SNPs_in_PAR]) #in genon finds the names of these markers
SNPs_to_remove<-append(SNPs_to_remove,markers_in_PAR)


#test for behaviour of SNPs on the X chromosome 
#males should not be heterozygous for SNPs on the X chromosome
#if a X SNP has many males appearing heterozygous for it, this SNP does not behave as we would expect 
#----------------------------------------------------------------

proportion_hetero_X<-apply(genon[male_rows,index_X_SNPs],2,function(x) sum(x=="1")/length(male_rows)) #for each SNP calculating proportion of heterozygous reads when only male samples
hetero_X_SNPs<-which(proportion_hetero_X>0.1) #from looking at graph of heterozgosity this looks to be a reasonable cut off but subsequent testing may mean its altered
markers_hetero_X<-colnames(genon[0,index_X_SNPs[hetero_X_SNPs]]) #finds the names of the markers to remove 
SNPs_to_remove<-append(SNPs_to_remove,markers_hetero_X) #adds name of markers to vector




#####one paper mentioned filtering via maf, "a lower maf means there is a smaller chance of seeing a heterozygous SNP and decreases the chance that a SNP will be identified and excluded" 
#####this paper got rid of SNPs with maf<0.015, verylow
#----------------------------------------------------------------

low_maf<-which(maf<0.015)
markers_low_maf<-colnames(genon[0,low_maf])
SNPs_to_remove<-append(SNPs_to_remove,markers_low_maf)

#check for missing data, check for SNPs on X chromosome that have more than 10% missing genotypes with female samples
#--------------------------------------------------------------
missing<-apply(genon[female_rows,index_X_SNPs],2,function(x) sum(is.na(x))/length(index_X_SNPs))
lot_missing<-which(missing>0.1)
rm_missing<-colnames(genon[0,index_X_SNPs[lot_missing]])

SNPs_to_remove<-append(SNPs_to_remove,rm_missing)

#once all checks are complete need to remove duplicates from the SNP list and make a csv file that can be loaded with subsequent datasets where sex is not verified
#----------------------------------------------------------------

SNPs_to_remove<-unique(SNPs_to_remove)
cat(length(SNPs_to_remove), "SNPs to remove\n") #prints out how many SNPs it will remove
cat(length(sex_markers$Marker.name)-length(SNPs_to_remove), "SNPs remain\n")

#code needed to create graphs
#----------------------------------------------------------------

cols_remove<-match(SNPs_to_remove,colnames(genon[0,])) #gets index of columns to remove from genon
modified_genon<-genon[,-cols_remove] #creates new genon without poorly behaving SNPs, prediction to be done from this matrix

modified_sex_markers<-sex_markers[-cols_remove,] #from file of sex markers remove the poorly behaving SNPs
modified_index_Y_SNPs<-which(modified_sex_markers$chromosome=="CM008042.1") #finds the new index of Y chromosome SNPs

cat(length(modified_index_Y_SNPs), "Y SNPs\n")
#modified_index_Y_SNPs<-modified_index_Y_SNPs[-11]
modified_index_X_SNPs<-which(modified_sex_markers$chromosome=="CM008041.1") #finds the new index of X chromosome SNPs

cat(length(modified_index_X_SNPs), "X SNPs\n")


#graph to look at heterozygoisty before PAR removed and before other filtering for males
#ideally all SNPs not in PAR region would show heterozygosity around zero
#----------------------------------------------------------------

SNP_pos_X<-sex_markers$SNP_position[index_X_SNPs] #creates vector with X SNP positions 
male_SNP_X_hetero<-apply(genon[male_rows,index_X_SNPs],2,function(x) sum(x=="1")/length(male_rows)) #calculates the heterozygosity the SNPs have with males, should be low
male_and_X<-data.frame(male_SNP_X_hetero,SNP_pos_X) #create dataframe 
male_and_X<-male_and_X[order(male_and_X$SNP_pos_X),] #orders the SNP positions smallest to largest 


png(file=file.path(output_file,"heterozygosity_on_X_males.png"), height=600,width=600)
plot(male_and_X$SNP_pos_X, male_and_X$male_SNP_X_hetero, xlab="X chr position", ylab="prprtn heterozygosity X chr", pch=19,cex.lab=1.5,cex.axis=1.5)
dev.off()

#graph to look at heterozygosity after PAR removed and filtering done for males
#after filtering we would expect this graph to show no SNPs with heterozygosity above 0.1
#ideally most SNPs will be at zero
#----------------------------------------------------------------

modified_SNP_pos_X<-modified_sex_markers$SNP_position[modified_index_X_SNPs]
modified_male_SNP_X_hetero<-apply(modified_genon[male_rows,modified_index_X_SNPs],2,function(x) sum(x=="1")/length(male_rows))
modified_male_and_X<-data.frame(modified_male_SNP_X_hetero,modified_SNP_pos_X)
modified_male_and_X<-modified_male_and_X[order(modified_male_and_X$modified_SNP_pos_X),] #orders the SNP positions smallest to largest 

png(file=file.path(output_file,"heterozygosity_on_X_males_modified.png"), height=600,width=600)
plot(modified_male_and_X$modified_SNP_pos_X, modified_male_and_X$modified_male_SNP_X_hetero, xlab="X chromosome position", ylab="proportion heterozygosity on X chromosome", ylim=c(0,0.7),main="Proportion heterozygoisty for SNPs on X chromosome and males after filtering",pch=19)
dev.off()

### For paper
# png("FigureS1.png",width=12*300,height=6*300, type='cairo', res=300)
# par(mfrow=c(1,2), mar=c(6.6,4.6,2,1))
# plot(male_and_X$SNP_pos_X/(1e6), male_and_X$male_SNP_X_hetero, xlab="X Position (Mb)",
#      ylab="Heterozygosity (X chr)", pch=19,cex.lab=1.5,cex.axis=1.5)
# abline(v=170, col=2, lwd=2)
# abline(h=0.1,col=2, lwd=2, lty=2)
# abline(h=0, lty=2)
# mtext("(a)", side=1, line=5,cex=1.5)
# plot(modified_male_and_X$modified_SNP_pos_X/(1e6), modified_male_and_X$modified_male_SNP_X_hetero, xlab="X Position (Mb)", ylab="Heterozygosity (X chr)", ylim=c(0,0.7175573),main="",pch=19,cex.lab=1.5,cex.axis=1.5)
# abline(h=0, lty=2)
# mtext("(b)", side=1, line=5,cex=1.5)
# dev.off()


#graph to look at heterozygoisty before PAR removed and other filtering for females
#would expect a fairly uniform distribution of data points
#----------------------------------------------------------------

SNP_pos_X_females<-sex_markers$SNP_position[index_X_SNPs] #creates vector with X SNP positions 
female_SNP_X_hetero<-apply(genon[female_rows,index_X_SNPs],2,function(x) sum(x=="1")/length(female_rows)) #calculates the heterozygosity the SNPs have with females
female_and_X<-data.frame(female_SNP_X_hetero,SNP_pos_X_females) #create dataframe 
female_and_X<-female_and_X[order(female_and_X$SNP_pos_X_females),] #orders the SNP positions smallest to largest 

png(file=file.path(output_file,"heterozygosity_on_X_females.png"), height=600,width=600)
plot(female_and_X$SNP_pos_X, female_and_X$female_SNP_X_hetero, xlab="X chromosome position", ylab="proportion heterozygosity on X chromosome", main="Proportion heterozygoisty for SNPs on X chromosome and females",pch=19)
dev.off()

#graph to look at heterozygosity after PAR removed and filtering done for females
#would expect a fairly similar distribution except SNPs in PAR are missing as well as SNPs
#filtered out for showing high amounts of heterozygosity in males 
#----------------------------------------------------------------

modified_SNP_pos_X_females<-modified_sex_markers$SNP_position[modified_index_X_SNPs]
modified_female_SNP_X_hetero<-apply(modified_genon[female_rows,modified_index_X_SNPs],2,function(x) sum(x=="1")/length(female_rows))
modified_female_and_X<-data.frame(modified_female_SNP_X_hetero,modified_SNP_pos_X_females)
modified_female_and_X<-modified_female_and_X[order(modified_female_and_X$modified_SNP_pos_X),] #orders the SNP positions smallest to largest 

png(file=file.path(output_file,"heterozygosity_on_X_females_modified.png"), height=600,width=600)
plot(modified_female_and_X$modified_SNP_pos_X_females, modified_female_and_X$modified_female_SNP_X_hetero, xlab="X chromosome position", ylab="proportion heterozygosity on X chromosome", ylim=c(0,0.7),main="Proportion heterozygoisty for SNPs on X chromosome and females after filtering",pch=19)
dev.off()

#graph to look at Y chromosome and male heterozygosity
#would expect all SNPs on Y for males to show low amounts of heterozygosity
#potentially if SNPs with high heterozygosity are found near the end of the chromosome
#it could indicate a SNP in the PAR region 
#----------------------------------------------------------------

SNP_pos_Y<-sex_markers$SNP_position[index_Y_SNPs]
male_SNP_Y_hetero<-apply(genon[male_rows,index_Y_SNPs],2,function(x) sum(x=="1")/length(male_rows))
male_and_Y<-data.frame(SNP_pos_Y,male_SNP_Y_hetero)
male_and_Y<-male_and_Y[order(male_and_Y$SNP_pos_Y),]

png(file=file.path(output_file,"heterozygosity_on_Y_males.png"),height=600,width=600)
plot(male_and_Y$SNP_pos_Y,male_and_Y$male_SNP_Y_hetero,xlab="Y chromosome position",ylab="proportion heterozygosity on Y",main="proportion heterozygosity for males on Y at each marker",pch=19)
dev.off() 

#graph to look at Y chromosome after filtering with heterozygosity 
#would expect this graph to have SNPs showing little to no heterozygosity
#----------------------------------------------------------------

modified_SNP_pos_Y<-modified_sex_markers$SNP_position[modified_index_Y_SNPs]
modified_male_SNP_Y_hetero<-apply(modified_genon[male_rows,modified_index_Y_SNPs],2,function(x) sum(x=="1")/length(male_rows))
modified_male_and_Y<-data.frame(modified_SNP_pos_Y,modified_male_SNP_Y_hetero)
modified_male_and_Y<-modified_male_and_Y[order(modified_male_and_Y$modified_SNP_pos_Y),]

png(file=file.path(output_file,"heterozygosity_on_Y_males_modified.png"),height=600,width=600)
plot(modified_male_and_Y$modified_SNP_pos_Y,modified_male_and_Y$modified_male_SNP_Y_hetero,xlab="Y chromosome position",ylab="proportion heterozygosity on Y",ylim=c(0,0.85),main="proportion heterozygosity for males on Y at each marker after filtering",pch=19)
dev.off() 

#graph to look at females and the SNPs on Y that record hits
#would expect that females have no readins for SNPs on the Y chromosome
#----------------------------------------------------------------

SNP_positions_Y<-sex_markers$SNP_position[index_Y_SNPs]
females_and_Y_SNPs<-apply(genon[female_rows,index_Y_SNPs],2, function(x) sum(!is.na(x))/length(female_rows)) #should be zero 
females_and_Y<-data.frame(SNP_positions_Y,females_and_Y_SNPs)
females_and_Y<-females_and_Y[order(females_and_Y$SNP_positions_Y),]


png(file=file.path(output_file,"hits_on_Y_females.png"),height=600,width=600)
plot(females_and_Y$SNP_positions_Y,females_and_Y$females_and_Y_SNPs, xlab="Y chromosome position", ylab="proportion of SNPs on Y", main="Females and proportion of SNPs on Y chromosome",pch=19)
dev.off()

#graph to look at females and SNPs on Y that record hits after filtering 
#this graph should show all SNPs around zero or with very low readings 
#----------------------------------------------------------------

modified_SNP_positions_Y<-modified_sex_markers$SNP_position[modified_index_Y_SNPs]
modified_females_and_Y_SNPs<-apply(modified_genon[female_rows,modified_index_Y_SNPs],2, function(x) sum(!is.na(x))/length(female_rows)) #should be zero 
modified_females_and_Y<-data.frame(modified_SNP_positions_Y,modified_females_and_Y_SNPs)
modified_females_and_Y<-modified_females_and_Y[order(modified_females_and_Y$modified_SNP_positions_Y),]

png(file=file.path(output_file,"hits_on_Y_females_modified.png"),height=600,width=600)
plot(modified_females_and_Y$modified_SNP_positions_Y,modified_females_and_Y$modified_females_and_Y_SNPs, xlab="Y chromosome position", ylab="proportion of SNPs on Y",ylim=c(0,1.0), main="Females and proportion of SNPs on Y chromosome after filtering",pch=19)
dev.off()

#code for filtering the deer tagdigger and deerAll20170301 file to ony sex SNPs that got through the filtering
#----------------------------------------------------------------

cat("Plots completed, subsetting deer count file to only include filterd sex SNPs\n")


sex_tagpair_number<-subset(sex_SNPs, select=c("tagpair","chrom","pos"))

subset_tagdigger<-tag_digger[tag_digger$UNEAK %in% sex_tagpair_number$tagpair,] #matches UNEAK tagpair to tagpair of sex chromosomes in sex_SNP_pos file

names<-SNPs_to_remove
split_names<-strsplit(as.character(names),"_")
marker_name<-lapply(split_names, "[[",1)
marker_names_to_remove<-unlist(marker_name, use.names=FALSE)
names_to_remove<-which(subset_tagdigger$Marker.name %in% marker_names_to_remove)
subset_tagdigger<-subset_tagdigger[-names_to_remove,]

sex_markers<-subset_tagdigger[,"Marker.name"] #makes vector of just the marker names


deer_count_header<-fread(deer_count,nrows=0,data.table=FALSE) #file name of the counts/reads in only the header row i.e the name of the markers 
marker_name<-colnames(deer_count_header) #creates a vector of marker names 
split_marker_name<-strsplit(as.character(marker_name), "_") #splits the name byunderscore
sex_marker_name<-lapply(split_marker_name, "[[",1) #takes the first part of the string split
vector_marker_name<-unlist(sex_marker_name, use.names=FALSE) #turns into a vector
vector_cols<-which(vector_marker_name %in% sex_markers) #find sex markers index
deer_count_names<-c("V1",marker_name[vector_cols])#finds names of sex_markers, include first column 


#finding names of the markers that KGD removed initially
#---------------------------------------------------------------
colname<-colnames(genon)
find<-setdiff(deer_count_names,colname)
dif<-grep("_0",find)
f<-find[dif]

f<-append(f,SNPs_to_remove)
write.csv(f,file=file.path(output_file, "markers_to_remove.csv"),row.names=FALSE)


subset_deer<-fread(deer_count,select= deer_count_names,  showProgress=TRUE ) #reads only the columns of the of the sample IDs and sex markers/name of file needs to be same as above 
write.csv(subset_deer, file=file.path(output_file, "deerAll20170301_counts_filtered_sexmarkers_only.csv"), row.names = FALSE) #writes out file which is used in the next script, name of file is name of outfile you want head

#write.csv(deer_count_names, file=file.path(output_file, "markers_to_keep.csv"),row.names=FALSE)


cat("Code complete, final output file deerAll20170301_counts_filtered_sexmarkers_only.csv made\n")  

rm("allelecounts","alleles","depth","depth.orig","fcolo","genon","HWdis",
   "l10LRT","maf","nind","nsnps","p","pg","RAcounts","sampdepth","sampdepth.max","sampdepth.med",
   "seqID","SNP_Names","snpdepth", envir=globalenv())

}
