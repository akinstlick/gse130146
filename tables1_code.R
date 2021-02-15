setwd("D:/SingleCellAnalysis/GSE130146/Results/")
#--------------------------------------
#this document finds marker genes from table S1 provided by the authors
topvals = 10


tableS1 <- read.table("D:/SingleCellAnalysis/GSE130146/Data/EB/TableS1_markGene.txt",
                      sep = "\t",skip=2,head=TRUE)

clusterid=unique(tableS1$cluster)
clusterid
#[1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17

n = length(clusterid)
n
#[1] 18

summary(tableS1$p_val_adj)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00e+00 0.00e+00 0.00e+00 4.03e-02 9.80e-06 1.00e+00

markerS1 <- list()

for (i in clusterid){
  temp = subset(tableS1,cluster==i)
  temp = subset(temp,p_val_adj < .05)
  temp2 = order(temp$avg_logFC,decreasing = TRUE)
  temp = temp[temp2,][1:topvals,]
  markerS1[[i+1]] = temp$gene.1
}

clusterid[6] = "endoderm"
clusterid[11] = "hemangiogenic"
clusterid[14] = "smooth muscle"
clusterid[15] = "cardiac"
clusterid[16] = "primordial germ"
clusterid[18] = "naive pluripotent"

names(markerS1) = clusterid

#extracted from published table s1 the author's biomarkers for 18 clusters

#now: manually replace numerical ids with biological names

#then save markerS1 somewhere else

save(markerS1, file="markerS1.RData", compress=TRUE)