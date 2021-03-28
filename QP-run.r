library(SignatureEstimation)
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
QP=decomposeQP(mut1, cancer_signatures)
write.table(QP,'Expo.out')
