# Author: Pierre Merckaert
# Contact: merckaert.pierre@gmail.com
# Date: 06/2018
# 
# Description : Create the columns for the input dataframe (specific order to make the predictions)
# input : dataframe of the 30mer sequence 
# output : dataframe of the 30mer sequence  and its features to extract (GC%, Melting Temp, indep order 1&2, dependant order 1&2,NGGN)

createColumns = function(input){
  nt = c('A','C','G','T')
  dint = c('AA','AC','AG','AT',
           'CA','CC','CG','CT',
           'GA','GC','GG','GT',
           'TA','TC','TG','TT')
  
  #Create ordered columns for 23mer, GC content and melting temperatures
  input$gRNA_23mer = 0
  input$GCcont_20mer = 0
  input$Tm0_6 = 0
  input$Tm7_14 = 0
  input$Tm15_19 = 0
  input$Tm30mer = 0
  
  #create columns on single nt amount (dependant order1)
  for (j in seq(length(nt))) {
    #amount
    nbColPos = paste('nb',nt[j],sep = '')
    input[[nbColPos]] = 0
  }
  
  #create columns on single dint amount (dependant order2)
  for (j in seq(length(dint))) {
    #amount
    nb_ColPos = paste('nb',dint[j],sep = '')
    input[[nb_ColPos]] = 0
  }
  
  #create columns on single nt positions (independant order1)
  for (j in seq(length(nt))) {
    #position
    for (loc in seq(-4,25)) {
      ntColPos = paste(nt[j],loc,sep = '')
      input[[ntColPos]] = 0
    }
  }
  
  #create columns on di-nt positions (independant order2)
  for (j in seq(length(dint))) {
    #position
    for (loc in seq(-4,24) ){
      nt_ColPos = paste(dint[j],loc,sep = '')
      input[[nt_ColPos]] = 0
    }
  }
  
  #create columns on di-nt flanking PAM 'GG'
  for (j in seq(length(dint))) {
    nggn_ColPos = paste(substr(dint[j], 1,1),paste(substr(dint[j], 2,2),'20', sep=''), sep= 'GG')
    input[[nggn_ColPos]] = 0
  }
  return(input)
}