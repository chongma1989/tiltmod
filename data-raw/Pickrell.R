# source("https://www.bioconductor.org/biocLite.R")
# biocLite("tweeDEseqCountData")
library("tweeDEseqCountData")
library("edgeR")
library("limma")

#==============================================#
# obtain and filter the RNA-seq count data     #
#==============================================#
data(pickrell1)
Counts=exprs(pickrell1.eset)
Gender=pickrell1.eset$gender
data("annotEnsembl63")
annot=annotEnsembl63[,c("Symbol","Chr")]
y=DGEList(counts=Counts, genes=annot[rownames(Counts),], group = Gender)

# keep genes with at least 1 count-per-million 
# (cpm) in at least 20 samples
isexpr=rowSums(cpm(y)>1)>=20
# keep only genes with defined annotation
hasannot= rowSums(is.na(y$genes))==0
y=y[isexpr & hasannot, , keep.lib.sizes=FALSE]
y=calcNormFactors(y) ## calculate normal factors
save(y,file='data/pickrell.rda',compress='xz')

