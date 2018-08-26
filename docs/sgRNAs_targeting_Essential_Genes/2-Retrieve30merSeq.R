#retrieve the 30mer sequence from the exome by searching for the 23mer on each strand

library(stringr)
input = read.csv('REP2_ML_Test_Set.csv', sep = ',', stringsAsFactors = F)
guide = read.table('dmel-all-exon-r6.09.fasta',sep = ';', stringsAsFactors = F)
guiderevcomp = read.table('dmel-all-exon-r6.09_REVCOMP.fasta',sep = ',', stringsAsFactors = F)

sgRNA = input$sgRNA23mer
DNAplus = guide$V1
DNAminus = guiderevcomp$V1


for (l in seq(length(sgRNA))){
  #search on the exome
  searchPlus = grep(paste("....",sgRNA[l], "...", sep=""), DNAplus)
  if(length(searchPlus != 0)){
    MER = str_extract(DNAplus[searchPlus[1]], paste("....",sgRNA[l], "...", sep=""))
  } else {
    #search on the reverse complement of the exome (complement strand)
      searchMinus = grep(paste("....",sgRNA[l], "...", sep=""), DNAminus)
      if(length(searchMinus != 0)){
      MER = str_extract(DNAminus[searchMinus[1]], paste("....",sgRNA[l], "...", sep=""))
      } else{
          MER = 'ERROR'
        }
  }
  input$sgRNA30mer[l] = MER
  print(l)
}

write.csv(input,'REP2_ML_Test_Set.csv',row.names = F)
