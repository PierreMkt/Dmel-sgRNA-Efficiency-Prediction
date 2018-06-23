
createColumns = function(input){
  #Create ordered columns 
  input$gRNA_23mer = 0
  input$GCcont_20mer = 0
  input$Tm0_6 = 0
  input$Tm7_14 = 0
  input$Tm15_19 = 0
  input$Tm30mer = 0
  
  #create columns on single nt amount
  for (j in seq(length(nt))) {
    #amount
    nbColPos = paste('nb',nt[j],sep = '')
    input[[nbColPos]] = 0
  }
  
  #create columns on single dint amount
  for (j in seq(length(dint))) {
    #amount
    nb_ColPos = paste('nb',dint[j],sep = '')
    input[[nb_ColPos]] = 0
  }
  
  #create columns on single nt positions
  for (j in seq(length(nt))) {
    #position
    for (loc in seq(-4,25)) {
      ntColPos = paste(nt[j],loc,sep = '')
      input[[ntColPos]] = 0
    }
  }
  
  #create columns on di-nt positions
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

#Tm calc function
calcTm = function(seq){
  nbA = str_count(seq,'A')
  nbC = str_count(seq,'C')
  nbT = str_count(seq,'T')
  nbG = str_count(seq,'G')
  
  if(str_length(seq) <= 13){
    Tm = (nbA+nbT)*2 + (nbG+nbC)*4 
  }
  else{
    Tm= 64.9 +41*(nbG+nbC-16.4)/(nbA+nbT+nbG+nbC)
  }
  return(Tm)
}
