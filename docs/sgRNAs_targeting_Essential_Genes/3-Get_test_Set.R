# sgRNA efficiency is biased towards 60/70%
# so we select a subset for the test set that is not so biased 
# by selecting an equal amount of sgRNA in each efficiency percentile

REP2 = read.csv("REP2count_Normalized.csv", header = T, sep = ',')
plot(REP2$Log2FC_Norm)

low001 = subset(REP2,(REP2$Log2FC_Norm <=0.1))
low0102 = subset(REP2,(REP2$Log2FC_Norm > 0.1) & (REP2$Log2FC_Norm <=0.2))
low0203 = subset(REP2,(REP2$Log2FC_Norm > 0.2) & (REP2$Log2FC_Norm <=0.3))
low0304 = subset(REP2,(REP2$Log2FC_Norm > 0.3) & (REP2$Log2FC_Norm <=0.4))
low0405 = subset(REP2,(REP2$Log2FC_Norm > 0.4) & (REP2$Log2FC_Norm <=0.5))
low0506 = subset(REP2,(REP2$Log2FC_Norm > 0.5) & (REP2$Log2FC_Norm <=0.6))
low0607 = subset(REP2,(REP2$Log2FC_Norm > 0.6) & (REP2$Log2FC_Norm <=0.7))
low0708 = subset(REP2,(REP2$Log2FC_Norm > 0.7) & (REP2$Log2FC_Norm <=0.8))
low0809 = subset(REP2,(REP2$Log2FC_Norm > 0.8) & (REP2$Log2FC_Norm <=0.9))
low091 = subset(REP2,(REP2$Log2FC_Norm > 0.9) & (REP2$Log2FC_Norm <=1.0))

test_df = rbind(low001,
            low0102,
            low0203,
            subset(low0304, sgRNA %in% sample(unique(low0304$sgRNA), 150)),
            subset(low0405, sgRNA %in% sample(unique(low0405$sgRNA), 150)),
            subset(low0506, sgRNA %in% sample(unique(low0506$sgRNA), 150)),
            subset(low0607, sgRNA %in% sample(unique(low0607$sgRNA), 150)),
            subset(low0708, sgRNA %in% sample(unique(low0708$sgRNA), 150)),
            low0809,
            low091)

plot(test_df$Log2FC_Norm)


sgRNA_seq = read.csv("CRISPR_KO_library_chip_order.csv", header = T, sep = ',');sgRNA_seq = sgRNA_seq[order(sgRNA_seq$DsgID),]
test_seq = sgRNA_seq[sgRNA_seq$DsgID %in% test_df$sgRNA,];test_seq = test_seq[order(test_seq$DsgID),]
test_df =test_df[order(test_df$sgRNA),]
test_df$sgRNA23mer = test_seq$gRNA1

write.csv(test_df[,c(1,7,6)],"REP2_ML_Test_Set.csv",row.names = F)
