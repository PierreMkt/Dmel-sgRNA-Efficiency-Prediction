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
library(stringr)

args = commandArgs(trailingOnly=F)

#retrieve the location of the file from the arguments
PATH = str_extract(substr(args[4],8,str_length(args[4])),'.*/')

if(is.na(PATH)){
  PATH = './'
}

if(!exists("preprocess", mode="function")) source(str_c(PATH,"FeatureExtraction_Preprocessing.R"))
if(!exists("createColumns", mode="function")) source(str_c(PATH,"OrderedColumns_Preprocessing.R"))
if(!exists("argsProcessing", mode="function")) source(str_c(PATH,"OrderedColumns_Preprocessing.R"))

l_results = argsProcessing(args)
input = l_results$df
Nmer = l_results$Nmer

print('Extracting Features from raw the sequence')

#create the columns of the dataframe and extract the features from the input sequence(s) and fill the dataframe columns 
input = createColumns(input,Nmer)
input = preprocess(input,Nmer)
output = "R_Featurized_sgRNA.csv"

#output the dataframe of the extracted features
write.csv(input,str_c(PATH,output),row.names = F)
