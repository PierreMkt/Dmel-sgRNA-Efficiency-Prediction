#!/usr/bin/env Rscript

# Author: Pierre Merckaert
# Contact: merckaert.pierre@gmail.com
# Date: 06/2018
# 
# Description: Extract the features required for efficiency prediction from 30mer sgRNA sequences
# Input: 30mer sgRNA or csv file listing multiple sgRNA in the first column
# Output: csv file of all the features extracted from the input sequences

#install required packages and import functions
if(!require(stringr)) install.packages("stringr",repos = "http://cran.us.r-project.org")
if(!exists("preprocess", mode="function")) source("Rpreprocessing/FeatureExtraction_Preprocessing.R")
if(!exists("createColumns", mode="function")) source("Rpreprocessing/OrderedColumns_Preprocessing.R")

args = commandArgs(trailingOnly=TRUE)
print('Extracting Features from input')
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input 30mer sequence OR csv file)", call.=FALSE)
} else if (length(args)==1) {
    #error if input file is not a .csv 
    if (substr(args[1],str_length(args[1])-3,str_length(args[1])) == '.csv'){
      input = read.csv(args[1], sep = ',', stringsAsFactors = F,header = T)
      colnames(input) = c("gRNA_30mer")
    }
    #create the dataframe if the input is a 30nt sequence
    else if(str_length(args[1])==30){
      input = data.frame(gRNA_30mer=args[1])
    }
  #error if input sequence is not 30 nucleotides long
    else{
      stop("input sequence must be 30 nucleotides long (-4 to 25, with NGG PAM at position 20)", call.=FALSE)
    }
}

#create the columns of the dataframe
input = createColumns(input)

#extract the features from the input sequence(s) and fill the dataframe columns 
input = preprocess(input)

#output the dataframe of the extracted features
write.csv(input,"R_Featurized_sgRNA.csv",row.names = F)
