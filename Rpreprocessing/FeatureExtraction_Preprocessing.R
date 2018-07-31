# Author: Pierre Merckaert
# Contact: merckaert.pierre@gmail.com
# Date: 06/2018
# 
# Description : Extract the features from the 30nt sequence(s) and fill the input dataframe
# input : dataframe of the 30mer sequence and its features to extract (GC%, Melting Temp, indep order 1&2, dependant order 1&2,NGGN)
# output : dataframe filled with the 30mer sequence and its extracted features
preprocess = function(input,Nmer){
  nt = c('A','C','G','T')
  dint = c('AA','AC','AG','AT',
           'CA','CC','CG','CT',
           'GA','GC','GG','GT',
           'TA','TC','TG','TT')
  
  if(Nmer==30){
    start = 5
    end = 24
  } else {
    start = 1
    end = 20
  }
  
  # Progression Bar setup
  total <- length(input[,1])
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  input[,1] = toupper(input[,1])
  
  for (l in seq(length(input[,1]))) {
    # Progression Bar increment
    if(l %% 100 == 0.0){
      setTxtProgressBar(pb, l)
    }
    
    #GC content
    nbG = str_count(substr(input[,1][l],start,end),'G')
    nbC = str_count(substr(input[,1][l],start,end),'C')
    GCcont = (nbG+nbC)/str_length(substr(input[,1][l],start,end))
    input$GCcont_20mer[l] = GCcont
    
    #Tm
    if(Nmer==30){
      input$gRNA_23mer[l] = substr(input[,1][l],start,end+3)
      input$Tm0_6[l] = calcTm(substr(input[,1][l], start, start+6))
      input$Tm7_14[l] = calcTm(substr(input[,1][l], start+7, start+14))
      input$Tm15_19[l] = calcTm(substr(input[,1][l], start+15, end))
      input$Tm30mer[l] = calcTm(input[,1][l])
    }else if(Nmer==23){
      input$Tm1_7[l] = calcTm(substr(input[,1][l], start, start+6))
      input$Tm8_15[l] = calcTm(substr(input[,1][l], start+7, start+14))
      input$Tm16_20[l] = calcTm(substr(input[,1][l], start+15, end))
      input$Tm23mer[l] = calcTm(input[,1][l])
    }else if(Nmer==20){
      input$Tm1_7[l] = calcTm(substr(input[,1][l], start, start+6))
      input$Tm8_15[l] = calcTm(substr(input[,1][l], start+7, start+14))
      input$Tm16_20[l] = calcTm(substr(input[,1][l], start+15, end))
      input$Tm20mer[l] = calcTm(input[,1][l])
    }
    
    #populate the single nt position and amount columns (Order1)
    for (j in seq(length(nt))) {
      pos = unlist(gregexpr(nt[j], substr(input[,1][l], 1, str_length(input[,1][l]))))
      if(pos[1] != -1){
        for (k in seq(length(pos))) {
          if(Nmer==30){
            nt_ColPos=paste(nt[j],toString(pos[k]-5),sep = '')
          } else {nt_ColPos=paste(nt[j],toString(pos[k]),sep = '')}
          input[[nt_ColPos]][l] = 1
        }
        nb_ColPos=paste('nb',nt[j],sep = '')
        input[[nb_ColPos]][l] = length(pos)
      }
    }
    #populate the di-nt position and amount  columns (Order2)
    for (j in seq(length(dint))) {
      #pattern = (?=di-nt) and perl=TRUE to locate joined matches (CCC = 2 matches for CC di-nt)
      pos = unlist(gregexpr(paste("(?=",dint[j],")",sep=''), substr(input[,1][l], 1, str_length(input[,1][l])), perl=TRUE))
      if(pos[1] != -1){
        for (k in seq(length(pos))) {
          if(Nmer==30){
            nt_ColPos=paste(dint[j],toString(pos[k]-5),sep = '')
          } else {nt_ColPos=paste(dint[j],toString(pos[k]),sep = '')}
          input[[nt_ColPos]][l] = 1
        }
        nb_ColPos=paste('nb',dint[j],sep = '')
        input[[nb_ColPos]][l] = length(pos)
      }
      
      #populate the di-nt NGGN at position 20 column
      if(Nmer==30){
        if( (input[[paste(substr(dint[j], 1,1),'20',sep='')]][l]==1) & (input[[paste(substr(dint[j], 2,2),'23',sep='')]][l]==1) ){
          nggnColPos = paste(substr(dint[j], 1,1),paste(substr(dint[j], 2,2),'20', sep=''), sep= 'GG')
          input[[nggnColPos]][l] = 1
        }
      }
    }
  }
  setTxtProgressBar(pb, l)
  close(pb)
  return(input)
}

#function to calculate the melting temperature
calcTm = function(seq){
  nbA = str_count(seq,'A')
  nbC = str_count(seq,'C')
  nbT = str_count(seq,'T')
  nbG = str_count(seq,'G')
  
  if(str_length(seq) <= 13){
    Tm = (nbA+nbT)*2 + (nbG+nbC)*4 
  }  else {
    Tm= 64.9 +41*(nbG+nbC-16.4)/(nbA+nbT+nbG+nbC)
  }
  return(Tm)
}
