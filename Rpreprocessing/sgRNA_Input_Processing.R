#!/usr/bin/env Rscript
if(!require(stringr)) install.packages("stringr",repos = "http://cran.us.r-project.org")
if(!exists("preprocess", mode="function")) source("Rpreprocessing/Preprocessing.R")
if(!exists("createColumns", mode="function")) source("Rpreprocessing/util.R")


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input 30mer sequence)", call.=FALSE)
} else if (length(args)==1) {
    # if(str_length(args[1])!=30){
    #   stop("input sequence must be 30 nucleotides long (-4 to 25, with NGG PAM at position 20)", call.=FALSE)
    # }
    # default output file
    args[2] = "Processed_Sequence.csv"
}

nt = c('A','C','G','T')
dint = c('AA','AC','AG','AT',
         'CA','CC','CG','CT',
         'GA','GC','GG','GT',
         'TA','TC','TG','TT')

input = read.csv(args[1], sep = ',', stringsAsFactors = F,header = T)
# input = data.frame(gRNA_30mer = args[1])

input = createColumns(input)

input = preprocess(input)

write.csv(input,args[2],row.names = F)
