REP2_reads = read.table("REP2count.csv", sep = ',', header = T)
REP2_reads$rep2.read.count..45d = as.integer(REP2_reads$rep2.read.count..45d);
REP2_reads$REF.read.count = as.integer(REP2_reads$REF.read.count);
#Remove the 10% lowest INITIAL read counts
REP2_reads = REP2_reads[order(REP2_reads$REF.read.count),]
tenth = as.integer(10*length(REP2_reads$Gene)/100)
REP2_reads = REP2_reads[tenth:length(REP2_reads$Gene),]

REP2_reads = REP2_reads[,c(1,2,3,4)]
REP2_norm = REP2_reads

#NORMALIZATION
#1)Log of read counts
REP2_norm$LOGREF.read.count = log(REP2_norm$REF.read.count)
REP2_norm$LOGrep2.read.count..45d = log(REP2_norm$rep2.read.count..45d)

#2)Average each row
REP2_norm$LogAvg = (REP2_norm$LOGREF.read.count + REP2_norm$LOGrep2.read.count..45d)/2

#3)Filter genes with AvgLog = 0
REP2_norm = REP2_norm[order(REP2_norm$LogAvg),]
REP2_norm = REP2_norm[!is.infinite(REP2_norm$LogAvg),]
REP2_norm = REP2_norm[!REP2_norm$LogAvg==0,]

#4)Substract the AvgLog from the logCounts
REP2_norm$LOGREF.read.count = REP2_norm$LOGREF.read.count - REP2_norm$LogAvg
REP2_norm$LOGrep2.read.count..45d = REP2_norm$LOGrep2.read.count..45d - REP2_norm$LogAvg

#5&6)Calculate the median for each sample and convert to 'normal' numbers
medianREF = exp(median(REP2_norm$LOGREF.read.count))
median45 = exp(median(REP2_norm$LOGrep2.read.count..45d))

#7)Divide the original read counts by the scaling factors
REP2_norm$Rcounts_NORM = REP2_norm$REF.read.count/medianREF
REP2_norm$Rcounts45_NORM = REP2_norm$rep2.read.count..45d/median45

#8) Compute log2FC and [0:1] normalization
REP2_norm$Log2FC = log2(REP2_norm$Rcounts45_NORM / REP2_norm$Rcounts_NORM)
plot(REP2_norm$Log2FC)


for( i in seq(length(REP2_norm$Log2FC))){
  print(i)
  REP2_norm$Log2FC_Norm[i] = 1-((REP2_norm$Log2FC[i]-min(REP2_norm$Log2FC))/(max(REP2_norm$Log2FC)-min(REP2_norm$Log2FC)))
}

plot(REP2_norm$Log2FC)

write.csv(REP2_norm[,c(1,2,8,9,10,11)],"REP2count_Normalized.csv",row.names = F)



###########################################################
#Get sgRNA targeting selected essential genes
all_REP2sgRNA = read.csv("REP2count_Normalized.csv", header = T, sep = ',')
essential = read.csv("essGenes_REP2.csv", header = T, sep = ',')
sgRNA_seq = read.csv("CRISPR_KO_library_chip_order.csv", header = T, sep = ',')

ess_sgRNA = all_REP2sgRNA[ all_REP2sgRNA$Gene %in% essential$Gene ,];ess_sgRNA = ess_sgRNA[order(ess_sgRNA$sgRNA),]

ess_seq = sgRNA_seq[sgRNA_seq$DsgID %in% ess_sgRNA$sgRNA,];ess_seq = ess_seq[order(ess_seq$DsgID),]
ess_sgRNA$sgRNA23mer = ess_seq$gRNA1
ess_sgRNA = ess_sgRNA[order(ess_sgRNA$sgRNA), c(1,2,7,6)]

write.csv(ess_sgRNA, "7514sgRNA_ess_REP2.csv",row.names = F)
