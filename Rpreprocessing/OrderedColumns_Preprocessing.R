# Author: Pierre Merckaert
# Contact: merckaert.pierre@gmail.com
# Date: 06/2018
# 


argsProcessing = function(args){
  #error if input file is not a .csv 
  if (substr(args[6],str_length(args[6])-3,str_length(args[6])) == '.csv'){
    input = read.csv(args[6], sep = ',', stringsAsFactors = F,header = T)
    #check if the first row is a 20mer, 30mer or 23mer
    if(str_length(input[1,])==30){
      Nmer = 30
    } else if(str_length(input[1,])==23){
      Nmer = 23
    } else if(str_length(input[1,])==20){
      Nmer = 20
    } else {
      stop("Error : sgRNA length should be 20, 23 or 30 nucleotides long (see --help) ")
    }
    #iterate over all rows to check if they all are 20, 23 or 30 mer 
    for (i in seq(length(input[,]))) {
      if(str_length(input[i,])!=Nmer){
        stop(str_c("Error : row ",i," sgRNA length should be 20, 23 or 30 nucleotides long (see --help) "))
      }
    }
    colnames(input) = str_c("gRNA_",Nmer,"mer")
    
  } else if(str_length(args[6])==30){  #create the dataframe if the input is a 30nt sequence
    input = data.frame(gRNA_30mer=args[6])
    Nmer = 30
  } else if(str_length(args[6])==23){  #create the dataframe if the input is a 30nt sequence
    input = data.frame(gRNA_23mer=args[6])
    Nmer = 23
  } else if(str_length(args[6])==20){  #create the dataframe if the input is a 30nt sequence
    input = data.frame(gRNA_20mer=args[6])
    Nmer = 20
  } else{ #error if input sequence is not 30 nucleotides long
    stop("input sequence must be 30 nucleotides long (-4 to 25, with NGG PAM at position 20) or 20 (no PAM) or 23 nucleotides long", call.=FALSE)
  }
  l_results = list("df" = input,"Nmer" = Nmer)
  return(l_results)
}


# Description : Create the columns for the input dataframe (specific order to make the predictions)
# input : dataframe of the 30mer sequence 
# output : dataframe of the 30mer sequence  and its features to extract (GC%, Melting Temp, indep order 1&2, dependant order 1&2,NGGN)
createColumns = function(input,Nmer){
  nt = c('A','C','G','T')
  dint = c('AA','AC','AG','AT',
           'CA','CC','CG','CT',
           'GA','GC','GG','GT',
           'TA','TC','TG','TT')
  
  #Create ordered columns for 23mer, GC content and melting temperatures
  if(Nmer==30){
    input$gRNA_23mer = 0
    input$GCcont_20mer = 0
    input$Tm0_6 = 0
    input$Tm7_14 = 0
    input$Tm15_19 = 0
    input$Tm30mer = 0
    
    start = -4
    end = 25
  }else if(Nmer==23){
    input$GCcont_20mer = 0
    input$Tm1_7 = 0
    input$Tm8_15 = 0
    input$Tm16_20 = 0
    input$Tm23mer = 0
    
    start = 1
    end = 23
  }else{
    input$GCcont_20mer = 0
    input$Tm1_7 = 0
    input$Tm8_15 = 0
    input$Tm16_20 = 0
    input$Tm20mer = 0
    
    start = 1
    end = 20
  }
  
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
    for (loc in seq(start,end)) {
      ntColPos = paste(nt[j],loc,sep = '')
      input[[ntColPos]] = 0
    }
  }
  
  #create columns on di-nt positions (independant order2)
  for (j in seq(length(dint))) {
    #position
    for (loc in seq(start,end-1) ){
      nt_ColPos = paste(dint[j],loc,sep = '')
      input[[nt_ColPos]] = 0
    }
  }
  
  if(Nmer==30){
    #create columns on di-nt flanking PAM 'GG'
    for (j in seq(length(dint))) {
      nggn_ColPos = paste(substr(dint[j], 1,1),paste(substr(dint[j], 2,2),'20', sep=''), sep= 'GG')
      input[[nggn_ColPos]] = 0
    }
  }
  return(input)
}