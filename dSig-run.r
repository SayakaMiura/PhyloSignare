library(deconstructSigs)
mut=read.delim('MUTCOUNTIN',sep=',',header=FALSE,check.names=FALSE) #C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\Example\\E.csv
colnames(mut) <- c('Type','Count')
Pro=mut$Count/sum(mut$Count)
mut0=cbind(mut,Pro)
row.names(mut0)=mut0$Type
signaturesCOSMIC=read.delim('SIGIN',sep=',')  #C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\PhyloSig1\\cosmicSig.csv
row.names(signaturesCOSMIC)=signaturesCOSMIC$X
drops <- c("X")
signaturesCOSMIC=signaturesCOSMIC[ , !(names(signaturesCOSMIC) %in% drops)]
order = match(row.names(mut0), row.names(signaturesCOSMIC))
signaturesCOSMIC=signaturesCOSMIC[,c('SIGLS')]
cancer_signatures = signaturesCOSMIC[order,]
cancer_signatures=data.matrix(cancer_signatures, rownames.force = NA)
mut1=data.frame(mut0$Pro)
MutProT=t(mut0$Pro)
colnames(MutProT) <-mut0$Type
rownames(MutProT) <- c('A')
B=data.frame(MutProT,check.names=FALSE)
cancer_signaturesT=t(cancer_signatures)
cancer_signaturesT=t(cancer_signatures)
SigForm=as.matrix(cancer_signaturesT)
SigForm1=data.frame(SigForm,check.names=FALSE)
test = whichSignatures(tumor.ref = B,sample.id = "A",signatures.ref = SigForm1, signature.cutoff = 0.001,contexts.needed = FALSE)
dSig=t(test$weights)
write.table(dSig,'Expo.out')

