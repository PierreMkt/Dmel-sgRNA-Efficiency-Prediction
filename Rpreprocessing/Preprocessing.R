preprocess = function(input){
  total <- length(input$gRNA_30mer)

  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (l in seq(length(input$gRNA_30mer))) {
    
    if(l %% 100 == 0.0){
      setTxtProgressBar(pb, l)
    }
    
    input$gRNA_23mer[l] = substr(input$gRNA_30mer[l],5,27)
    
    # if(input$sgRNA_strand[l] == 1){
    #   input$gRNA_30mer[l] = str_extract(input$DNAplus[l], paste("....",input$gRNA1[l], "...", sep=""))
    #   } else{
    #     input$gRNA_30mer[l] = str_extract(input$DNAminus[l], paste("....",input$gRNA1[l], "...", sep=""))
    #   }

    #GC content
    nbG = str_count(substr(input$gRNA_30mer[l],5,24),'G')
    nbC = str_count(substr(input$gRNA_30mer[l],5,24),'C')
    GCcont = (nbG+nbC)/str_length(substr(input$gRNA_30mer[l],5,24))
    input$GCcont_20mer[l] = GCcont
    
    #Tm
    input$Tm0_6[l] = calcTm(substr(input$gRNA_30mer[l], 5, 11))
    input$Tm7_14[l] = calcTm(substr(input$gRNA_30mer[l], 12, 19))
    input$Tm15_19[l] = calcTm(substr(input$gRNA_30mer[l], 20, 24))
    input$Tm30mer[l] = calcTm(substr(input$gRNA_30mer[l], 1, str_length(input$gRNA_30mer[l])))
    
    #populate the single nt position and amount  columns
    for (j in seq(length(nt))) {
      pos = unlist(gregexpr(nt[j], substr(input$gRNA_30mer[l], 1, str_length(input$gRNA_30mer[l]))))
      if(pos[1] != -1){
        for (k in seq(length(pos))) {
          nt_ColPos=paste(nt[j],toString(pos[k]-5),sep = '')
          input[[nt_ColPos]][l] = 1
        }
        nb_ColPos=paste('nb',nt[j],sep = '')
        input[[nb_ColPos]][l] = length(pos)
      }
    }
    #populate the di-nt position and amount  columns
    for (j in seq(length(dint))) {
      
      #pattern = (?=di-nt) and perl=TRUE to locate joint matches (CCC = 2 matches for CC di-nt)
      pos = unlist(gregexpr(paste("(?=",dint[j],")",sep=''), substr(input$gRNA_30mer[l], 1, str_length(input$gRNA_30mer[l])), perl=TRUE))
      if(pos[1] != -1){
        for (k in seq(length(pos))) {
          nt_ColPos=paste(dint[j],toString(pos[k]-5),sep = '')
          input[[nt_ColPos]][l] = 1
        }
        nb_ColPos=paste('nb',dint[j],sep = '')
        input[[nb_ColPos]][l] = length(pos)
      }
      
      #populate the di-nt NGGN at position 20 column
      if( (input[[paste(substr(dint[j], 1,1),'20',sep='')]][l]==1) & (input[[paste(substr(dint[j], 2,2),'23',sep='')]][l]==1) ){
        nggnColPos = paste(substr(dint[j], 1,1),paste(substr(dint[j], 2,2),'20', sep=''), sep= 'GG')
        input[[nggnColPos]][l] = 1
      }
    }
  }
  setTxtProgressBar(pb, l)
  close(pb)
  return(input)
}