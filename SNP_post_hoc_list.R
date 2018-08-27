load("PAWS_methlylumi.RData")

SNPCpG<-fData(PAWS.2)[,c(5,8)]
save(SNPCpG, file="SNPCpG.RData")