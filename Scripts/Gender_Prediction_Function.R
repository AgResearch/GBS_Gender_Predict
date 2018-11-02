

#this function will predict the gender of samples loaded
#by running them and the selected sex markers through KGD
#will produce output file of samples with predicted and
#phenotyped genders
#also produces a file with samples that have different recorded
#and predicted genders
#--------------------------------------------------------------

predict_gender<-function(sampleID,deer_results,sex_markers,filtered_markers, output_file="./", withPlotly=F, group=NULL){

  if(!require(ggplot2))
    stop("ggplot2 library is not available")
  if(!require(data.table))
    stop("data.table library is not available")

dir.create(output_file)
#output_file<-"/dataset/Predicting_Gender_Deer/active/Code/files_input_output_for_prediction"

samples<-read.csv(file=sampleID, header=TRUE)
vector_of_samples<-samples[,"SampleID"]
cat("Number of samples:", length(vector_of_samples))

#code to ensure only SNPs that made it through the original
#filtering step are included in subsequent analysis
#sometimes KGD will re-introduce previously filtered SNPs
#------------------------------------------------------------

markers_filtered<-read.csv(file=filtered_markers,header=TRUE)
markers_to_remove<-as.character(markers_filtered[,"x"])
marker<-strsplit(as.character(markers_to_remove), "_")
m_name<-lapply(marker, "[[",1)
vector_m_name<-unlist(m_name, use.names=FALSE)

deer_count<-fread(file=deer_results,data.table=FALSE, showProgress=TRUE)

colname<-colnames(deer_count)
split_name<-strsplit(as.character(colname), "_")
name<-lapply(split_name, "[[",1)
vector_name<-unlist(name, use.names=FALSE)

match<-which(vector_name %in% vector_m_name)
colname_remove<-colname[match]

deer_count<-deer_count[,-match]

SampleID<-deer_count$V1 # making vector of the SampleIDs of this #file
split_sample_ID<-strsplit(as.character(SampleID), "_") #breaks #sampleID up by splitting at the underscores
Sample_number<- lapply(split_sample_ID, "[[",1) #takes the first #part of the entire sampleID
vector_sampleID<-unlist(Sample_number, use.names=FALSE) #turns #this above list into a vector
vector_rows<-which(vector_sampleID  %in% vector_of_samples) #finds the sampleIDs in the large file by matching with the #sampleIDs found from script 03
deer_count_sub<-deer_count[vector_rows,1:ncol(deer_count)]
write.csv(deer_count_sub, file=file.path(output_file,"deer_count_subsetted_for_prediction.csv"), row.names = F) #only #taking #rows of the #samples we want in count file

cat(length(vector_rows), "samples found\n") #prints number of samples found, which may be greater than the number of individuals

#code to load in and run KGD
#----------------------------------------------------------------

gform <- "TagDigger"  ####Used for RA files

#gform <- "uneak"   ####Used for HapMap files

genofile<-file.path(output_file,"deer_count_subsetted_for_prediction.csv")
###Add location to Ra or HapMap file
sampdepth.thresh <- 0.3

source("../Scripts/GBS-Chip-Gmatrix.R",local=TRUE)

cat("KGD complete\n")


#code that will act to merge samples from the same individuals
#this code was taken from Genomnz and the full code is with Ken
#Dodds
#--------------------------------------------------------------


ped0<-read.csv(file=sampleID, header=TRUE) #need a file #that has information on animalID
SampleID <- read.table(text=seqID,sep="_",fill=TRUE,stringsAsFactors=FALSE)[,1]
pedpos <- match(SampleID,ped0$SampleID)
indsubset <- which(!is.na(pedpos))
animalID <- ped0$AnimalID[pedpos[indsubset]] #the variable after #$ needs to be changed to the column name that holds ID

system.time(mg <- mergeSamples(animalID) )
if(!is.null(group) && is.vector(group) && length(group)==1 && (group %in% colnames(ped0))){
  plotlygroup <- as.character(ped0[[group]][pedpos[indsubset]])
  plotlygroup <- plotlygroup[match(mg$mergeIDs,animalID)]
}
else{
  plotlygroup <- NULL
  group=NULL
}

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
cat("Sample merge complete\n")
cat( length(animalID), "animals after merge\n")

#Assigning markers in genon as on X or Y chromosome
#-----------------------------------------------------------

sex_marker_file<-read.csv(file=sex_markers, header=TRUE) #location of file with markers and/ coresponding chromosome ,made #in 02 script
markers<-match(SNP_Names, sex_marker_file$Marker.name) #creates a #vector of row indices/of the markers we want out of the #sex_marker file
sex_markers<-sex_marker_file[markers,] #holds only the marker names on the X and Y chromosome


index_Y_SNPs<-which(sex_markers$chromosome == "CM008042.1") #holds indicies of columns in the genon matrix which have the Y SNPs
index_X_SNPs<-which(sex_markers$chromosome == "CM008041.1")#holds #indicies of columns in genon matrix with the X SNPs

cat(length(SNP_Names),"SNPs After filtering\n")
cat(length(index_Y_SNPs), "Y chromosome SNPs\n")
cat(length(index_X_SNPs), "X chromosome SNPs\n")


#works out proporton of SNPs an individual has on the Y chromosome
#a male should score highly and a female should score lowly
#a female should have no SNPs on Y chromosome so a low score, #ideally zero
#----------------------------------------------------------------


proportion_SNPs_Y<-apply(genon[,index_Y_SNPs],1, function(x) sum(!is.na(x))/length(index_Y_SNPs)) #sums number non-missing SNPs #over total SNPs on Y
maxy <- length(index_Y_SNPs)

#working out proportion of heterozygosity on the X chromosome
#males should have a low score and females a relatively higher score
#this code counts the amount of 1's an individual has for SNPs on the X chromosome
#it then works out the proportion of heterozygosity by divding
#the total number of non-missing points a sample has for the SNPs #on X chromosome
#----------------------------------------------------------------


proportion_Hetero_X<-apply(genon[,index_X_SNPs],1, function(x) sum(x=="1")/sum(!is.na(x))) #######the genon matrix may have to be changed to hetero_X_genon

#code that adjusts heterozygosity for depth of the samples
#-----------------------------------------------------------

depth_X_matrix<-depth.orig[,index_X_SNPs] #cuts to matrix with depth values for SNPs on X chromosome
K_matrix<-(1-2*(0.5^depth_X_matrix))  #makes matrix of all the K values for the depths
K_matrix[K_matrix == -1]<- NaN    #turns values that are -1 (so depth of zero) to NaN
row_sum_K_matrix<-rowSums(K_matrix,na.rm=TRUE) #sums up the rows
total_points_X<-apply(genon[,index_X_SNPs],1,function(x) sum(!is.na(x))) #works out the number of non-missing points on the X chromosome
E<-row_sum_K_matrix/total_points_X  #divides total K of a row by the number of SNPs on that chromosome
new_prop_X<-proportion_Hetero_X/E #new depth adjusted proportion of heterozygosity used in subsequent analysis


#making a data.frame with #animalID,proportion_SNPs_Y,proportion_Hetero_X,sampdepth and an #empty gender column
#----------------------------------------------------------------


gender<- NA #need a place in the data frame for gender to go, at #the moment it's empty
prediction_data<-data.frame(animalID,proportion_SNPs_Y,new_prop_X,sampdepth,gender)#creates data.frame with data needed to predict gender with empty #gender column

#prediction algorithm
#will give an output file with animalID and the predicted gender
#----------------------------------------------------------------

predict_sex<-function(prop_Y,prop_Hetero_X,gender){
if(prop_Y > 20*((prop_Hetero_X)^2)+0.2){
gender<-"M"
return(gender)
}
else if(prop_Y<0.1+prop_Hetero_X){
gender<-"F"
return(gender)
}
#else if(prop_Y>0.2 && prop_Hetero_X<0.05){  #once you correct #for #sampdepth you may find that you don't need this condition, #tend #to be samples with low depth here
#gender<-"M"
#return(gender)
#}
else{
gender<-"Uncertain"
return(gender)
}
}


gender_prediction<-apply(prediction_data,1,function(x)  predict_sex(as.numeric(x[2]),as.numeric(x[3]),x[4])) #calls the function #to predict gender
females<-sum(gender_prediction=="F") #counts number of females #found
males<-sum(gender_prediction == "M") #counts numbers of males #found
uncertains<-sum(gender_prediction == "Uncertain") #counts number #of uncertains found
cat(females, "females predicted\n") #prints number of females #predicted
cat(males, "males predicted\n") #prints number of males predicted
cat(uncertains, "not able to be predicted\n") #prints number of #uncertains


#this part of code can compare the predicted gender with the phenotyped gender if there is a phenotype for the individuals
#----------------------------------------------------------------


samples<-read.csv(file=sampleID, header=TRUE)
animals<-match(prediction_data$animalID, samples$AnimalID)
phenotyped_gender<-samples[animals, "Sex"]
gender_output<-data.frame(animalID, gender_prediction,phenotyped_gender, new_prop_X,proportion_SNPs_Y,sampdepth)

if(!is.null(group) && is.vector(group) && length(group)==1 && (group %in% colnames(ped0))){
  gender_output <- cbind(gender_output,group=plotlygroup)
  names(gender_output)[ncol(gender_output)] <- group
}

write.csv(gender_output,file=file.path(output_file, "gender_prediction.csv"),row.names=FALSE)
#----------------------------------------------------------------


#code to find if any of the predicted genders do not match the #recorded genders
#----------------------------------------------------------------

difference_search<-function(animalID,predicted,phenotyped){ #a #function that looks for differences in the predicted and #phenotyped gender
if(!predicted == phenotyped){    #in the data.frame gender_ouput #checks the two coloumns
return("different")
}
}

look_for_differences<-apply(gender_output,1,function(x) difference_search(x[1],x[2],x[3])) #calls the function above and #holds the animalIDs that have differnces
differences<-which(look_for_differences == "different") #gets the #indicies of rows that have recorded differences in gender
cat( length(differences), "samples have different predicted and phenotyped genders\n")
query_gender<-gender_output[differences,]  #creates a new #data.frame that holds the individuals with differences in gender

write.csv(query_gender,file=file.path(output_file, "gender_queries.csv"),row.names=FALSE) #writes out the above data frame into a csv file



#for all graphs need a gender factor and need the number of individuals
#----------------------------------------------------------------

gender<-gender_output[,"gender_prediction"] #create a factor of gender to colour the points with
num_individuals<-length(animalID)
individuals<-seq(1, num_individuals, by=1)


#code to create plot of individuals and amount  of SNPs they have on the Y chromosome
#males should score highly and be in the upper group
#----------------------------------------------------------------

#recorded_index <- (gender_output$phenotyped_gender == "F")*19 + (gender_output$phenotyped_gender == "M")*17 + (!(gender_output$phenotyped_gender %in% c("F","M")))*15
newColor <- c("red","blue","darkgreen")[as.numeric(gender_output$phenotyped_gender)]

ds <- data.frame(IND=individuals, Ysnp = proportion_SNPs_Y*maxy, Pheno = gender_output$phenotyped_gender)

ggplot(ds, aes(x=IND, y=Ysnp, colour = Pheno)) + geom_point() + labs(y = "Y chr SNPs", x = 'Individuals') +
  scale_color_manual(values=c("red","blue","gray"))

ggsave(paste0(output_file,"/Individuals_and_SNPs_on_Y.png"),width=6, height=6)


#code to create plot showing the proportion of heterozyosity on the X chromosome for individuals
#----------------------------------------------------------------

#png(file=file.path(output_file,"Individuals_and_heterozygosity_on_X.png"),height=600,width=600)

ds <- data.frame(IND=individuals, HeteroX = new_prop_X, Pheno = gender_output$phenotyped_gender)

ggplot(ds, aes(x=IND, y=HeteroX, colour = Pheno)) + geom_point() + labs(y = "Heterozygosity (X chr)", x = 'Individuals') +
  scale_color_manual(values=c("red","blue","gray"))

ggsave(paste0(output_file,"/Individuals_and_heterozygosity_on_X.png"), height=6, width=6)

#dev.off()

#code that merges the two graphs above, should have males in the top left
#females should should be spread along the bottom, towards the right
#----------------------------------------------------------------

#png(file=file.path(output_file, "heterozygosity_vs_SNPs_on_Y.png"),height=600,width=600)

ds <- data.frame(Ysnp=proportion_SNPs_Y*maxy, HeteroX = new_prop_X, Pheno = gender_output$phenotyped_gender)
wrong <- as.character(gender_output$gender_prediction) != as.character(gender_output$phenotyped_gender)
wrong_ds <- ds[which(wrong),]
ds <- ds[which(!wrong),]

maxx <- max(0.35,gender_output$new_prop_X)
upperlimit <- function(x){ 20*x^2+0.2}
poly1 <- data.frame(x1 = c(0,maxx,maxx,0,0), y1=c(0.1,0.1+maxx,0,0,0.1)*maxy)
poly2 <- data.frame(x2 = c(seq(0,0.2,0.01),0,0), y2=c(upperlimit(seq(0,0.2,0.01)),1,0)*maxy)

ggplot(data=ds, aes(x=HeteroX, y=Ysnp, colour = Pheno)) + geom_point() + labs(x = "Heterozygosity (X chr)", y = "Y chr SNPs") +
  scale_color_manual(values=c("red","blue","gray")) + theme_bw() +
  geom_point(data = wrong_ds) +
  geom_polygon(data=poly2, aes(x=x2, y=y2), fill = rgb(0,0,1,alpha=0.1), show.legend = F, inherit.aes = F) +
  geom_polygon(data=poly1, aes(x=x1, y=y1), fill = rgb(1,0,0,alpha=0.1), show.legend = F, inherit.aes = F)

ggsave(paste0(output_file,"/heterozygosity_vs_SNPs_on_Y.png"), width=6, height=6)

#dev.off()

if(withPlotly){
  if(require(plotly)){

    gender_plotly <- as.character(gender)
    gender_plotly[which(gender_plotly=="F")] <- "Females"
    gender_plotly[which(gender_plotly=="M")] <- "Males"
    gender_plotly <- paste0("Predicted: ",gender_plotly)
    pheno_plotly <- as.character(gender_output$phenotyped_gender)
    pheno_plotly[which(pheno_plotly=="F")] <- "Females"
    pheno_plotly[which(pheno_plotly=="M")] <- "Males"
    pheno_plotly <- paste0("Phenotype: ",pheno_plotly)
    samp.info <- list(ID=gender_output$animalID, Phenotype=gender_output$phenotyped_gender, Predicted=gender_output$gender_prediction,
                      "Sample Depth"=format(round(gender_output$sampdepth,2), nsmall = 2))
    if(!is.null(group) && is.vector(group) && length(group)==1 && (group %in% colnames(ped0))){
      samp.info <- c(samp.info, list(group=as.character(gender_output[[group]])))
      names(samp.info)[length(samp.info)] <- group
    }
    hovertext <- apply(sapply(1:length(samp.info),function(x) paste0(names(samp.info)[x],": ",samp.info[[x]]),simplify = TRUE),1,paste0,collapse="<br>")
    curveFun <- function(x) 20*x^2 + 0.2
    xline1 <- c(seq(0,0.2,length.out=800),0,0); yline1 = c(curveFun(seq(0,0.2,length.out=800)),1,0)
    xulim <- min(max(0.3,ceiling(max(gender_output$new_prop_X)*100)/100) + 0.02,1)
    maxy <- max(gender_output$proportion_SNPs_Y)

    if(is.null(group)) {
      temp_p <- plot_ly(x=xline1, y=yline1, alpha=0.1, opacity=0.1,colors=c("red","blue")) %>%
        add_polygons(hoverinfo = "none", color = I("blue"),showlegend=F,hoveron="points") %>%
        add_polygons(x=c(0,xulim,xulim,0,0), y=c(0.1,0.1+xulim,0,0,0.1),
                     hoverinfo = "none", color = I("red"),showlegend=F,hoveron="points") %>%

        add_markers(y=gender_output$proportion_SNPs_Y*maxy,x=gender_output$new_prop_X, opacity=1, alpha=1,
                    color=as.factor(pheno_plotly), marker=list(size=8),
                    hoverinfo="text", text=apply(sapply(1:length(samp.info),function(x) paste0(names(samp.info)[x],": ",samp.info[[x]]),simplify = TRUE),1,paste0,collapse="<br>"))   %>%
        layout( xaxis=list(title="Heterozygosity (X chr)"), yaxis=list(title="Y chr SNPs"))
    } else{
      temp_p <- plot_ly(x=xline1, y=yline1, alpha=0.1, opacity=0.1,colors=c("red","blue")) %>%
        add_polygons(hoverinfo = "none", color = I("blue"),showlegend=F,hoveron="points") %>%
        add_polygons(x=c(0,xulim,xulim,0,0), y=c(0.1,0.1+xulim,0,0,0.1),
                     hoverinfo = "none", color = I("red"),showlegend=F,hoveron="points") %>%
        add_markers(y=gender_output$proportion_SNPs_Y*maxy,x=gender_output$new_prop_X, opacity=1, alpha=1,
                    color=as.factor(pheno_plotly), marker=list(size=8),
                    symbol=as.factor(paste0(group," :",gender_output[[group]])),
                    hoverinfo="text", text=apply(sapply(1:length(samp.info),function(x) paste0(names(samp.info)[x],": ",samp.info[[x]]),simplify = TRUE),1,paste0,collapse="<br>"))   %>%
        layout( xaxis=list(title="Heterozygosity (X chr)"), yaxis=list(title="Y chr SNPs"))
    }
    currentWD <- getwd()
    setwd(file.path(output_file))
    htmlwidgets::saveWidget(temp_p, "gender-plot.html")
    setwd(currentWD)

  }
}

cat("Plots complete\n")
cat("Code complete\n")

rm("allelecounts","alleles","depth","depth.orig","fcolo","genon","HWdis",
   "l10LRT","maf","nind","nsnps","p","pg","RAcounts","sampdepth","sampdepth.max","sampdepth.med",
   "seqID","SNP_Names","snpdepth", envir=globalenv())
}






