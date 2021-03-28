library(MutationalPatterns)
#filename <- "C:\\Users\\tuf78332\\Desktop\\CancerSoftware\\sigLASSO-master\\sigLASSO-master\\data\\signatures_probabilities.txt"
#'C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\PhyloSig1\\cosmicSig.csv'
#cancer_signatures <- read.table(filename, sep = "\t", header = TRUE)

#mut_mat=read.delim('MUTCOUNTIN',sep=',',header=FALSE,check.names=FALSE) #C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\Example\\E.csv
#colnames(mut_mat) <- c('Type','Count')
#order = match(mut_mat$Type, cancer_signatures$Somatic.Mutation.Type)
#cancer_signatures = cancer_signatures[order,]
#cancer_signatures = as.matrix(cancer_signatures[,4:33])
#colnames(cancer_signatures) = as.character(1:30)
#mut_mat1=data.frame(mut_mat$Count)
#fit_res <- fit_to_signatures(mut_mat1, cancer_signatures)
#write.table(fit_res$contribution/sum(mut_mat$Count),'Expo.out')


########
filename <-'SIGIN'#'C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\PhyloSig1\\cosmicSig.csv'
cancer_signatures <- read.table(filename, sep=',',header = TRUE)
cancer_signatures=cancer_signatures[,c('X','SIGLS')]
colnames(cancer_signatures)[colnames(cancer_signatures)=="X"] <- "Somatic.Mutation.Type"
mut_mat=read.delim('MUTCOUNTIN',sep=',',header=FALSE,check.names=FALSE) #C:\\Users\\tuf78332\\Desktop\\mutational-signature\\PhyloSig\\Example\\E.csv
colnames(mut_mat) <- c('Type','Count')
order = match(mut_mat$Type, cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[order,]
cancer_signatures = as.matrix(cancer_signatures[,2:ncol(cancer_signatures)])
colnames(cancer_signatures) = as.character(1:ncol(cancer_signatures))
mut_mat1=data.frame(mut_mat$Count)
fit_res <- fit_to_signatures(mut_mat1, cancer_signatures)
write.table(fit_res$contribution/sum(mut_mat$Count),'Expo.out')
