#!/usr/bin/env Rscript

# Author: Pierre Merckaert
# Contact: merckaert.pierre@gmail.com
# Date: 06/2018
# 
# Description: Extract the features required for efficiency prediction from 30mer sgRNA sequences
# Input: 30mer sgRNA or csv file listing multiple sgRNA in the first column
# Output: csv file of all the features extracted from the input sequences
# The format of the output csv is (number of columns) : 30mer(1), 23mer(1), GC%(1), Tm(4), indep Order1(4), indep Order2(16), dep Order1(120), dep Order2(464), NGGN(16)

#install required packages and import functions
if(!require(stringr)) install.packages("stringr",repos = "http://cran.us.r-project.org")

args = commandArgs(trailingOnly=F)
PATH = str_extract(args[4],'/.*/')

if(!exists("preprocess", mode="function")) source(str_c(PATH,"FeatureExtraction_Preprocessing.R"))
if(!exists("createColumns", mode="function")) source(str_c(PATH,"OrderedColumns_Preprocessing.R"))


#error if input file is not a .csv 
if (substr(args[6],str_length(args[6])-3,str_length(args[6])) == '.csv'){
  input = read.csv(args[6], sep = ',', stringsAsFactors = F,header = T)
  colnames(input) = c("gRNA_30mer")
} else if(str_length(args[6])==30){  #create the dataframe if the input is a 30nt sequence
  input = data.frame(gRNA_30mer=args[6])
}else{ #error if input sequence is not 30 nucleotides long
  stop("input sequence must be 30 nucleotides long (-4 to 25, with NGG PAM at position 20)", call.=FALSE)
}

print('Extracting Features from input')
#create the columns of the dataframe
input = createColumns(input)

#extract the features from the input sequence(s) and fill the dataframe columns 
input = preprocess(input)

#output the dataframe of the extracted features
write.csv(input,str_c(PATH,"R_Featurized_sgRNA.csv"),row.names = F)
