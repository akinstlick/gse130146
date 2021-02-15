#5/10
#trying to redo the whole process, using chapter 34 as an example.

#data loading

library(DropletUtils)
fnameeb <- file.path("D:/SingleCellAnalysis/GSE130146/Data/EB")
eb.sce <- read10xCounts(fnameeb, col.names=TRUE)

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(eb.sce), 
               keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(eb.sce) <- anno[match(rownames(eb.sce), anno$GENEID),]

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
chr.loc <- mapIds(ens.mm.v97, keys=rownames(eb.sce),
                  keytype="GENEID", column="SEQNAME")
is.mito <- which(chr.loc=="MT")
#snapshot date: 2020-4-27

#drop where total values are 0 in df
library(scater)
df <- perCellQCMetrics(eb.sce, subsets=list(Mito=is.mito))
df

save(eb.sce, file="Results/eb.sce.RData", compress=TRUE)

#these steps worked, produced an annotated df with mito percents done
#now i need to drop where mito% > 5 or fewer than 5000 genes


#------------------------------------------------------
#filtering

#for filtering the cells, we followed the author's methods.

qc.nexprs <- df$detected < 2000
qc.mito <- df$subsets_Mito_percent > 5
discard <- qc.nexprs | qc.mito
DataFrame(NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), Total=sum(discard))


#[1] 737280     10
length(qc.nexprs)
#[1] 737280
table(qc.nexprs)
qc.nexprs
#FALSE   TRUE 
#2202 735078 
table(qc.mito)
qc.mito
#FALSE   TRUE 
#355411  44246 
table(discard)
discard
#FALSE   TRUE 
#1731 735549 
eb.sce
#class: SingleCellExperiment 
#dim: 33456 737280 
#metadata(1): Samples
#assays(1): counts
#rownames(33456): ENSMUSG00000064842 ENSMUSG00000051951 ... ENSMUSG00000096730 ENSMUSG00000095742
#rowData names(2): ID Symbol
#colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ... TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
#colData names(2): Sample Barcode
#reducedDimNames(0):
#altExpNames(0):

filtered <- eb.sce[,!discard]

dim(filtered)
#[1] 33456  1731

save(filtered,file="Data/EB/ebfiltered.RData")

#mitochondrial control
#follow chapter 6 example, then apply
#more than 5% mitochondria reads or fewer than 2000 unique genes were filtered out - original authors
#do they have the batch effect correctly from paper
#do it 3 times and combine data together?
#correct for batch effect if they did it
#i need to annotate dataset correctly


#-----------------------------------
#normalization


#the author did not provide details for their normalization, so we used normalization by deconvolution
# and normalization by library size
library(scater)
libsf.filtered <- librarySizeFactors(filtered)
summary(libsf.filtered)

hist(log10(libsf.filtered), breaks=100, xlab="Log10[Size factor]", col='grey80')

library(scran)
set.seed(100)
clust.filtered <- quickCluster(filtered) 
table(clust.filtered)

deconv.filtered <- calculateSumFactors(filtered, cluster=clust.filtered)
summary(deconv.filtered)

plot(libsf.filtered, deconv.filtered, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16
     #col=as.integer(factor(sce.zeisel$level1class))
     #don't know classification, so we can't do this
                           )
abline(a=0, b=1, col="red")

head(rowData(filtered))

#----------------------------------------------------
#finding housekeeping genes

temp = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Mouse_atlas.txt", sep="\t",head=T)
temp2 = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Mouse_development.txt", sep="\t",head=T)
temp3 = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Human_development.txt", sep="\t",head=T)

temp <- subset(temp, Stability.index>0.8)$GeneSymbol
temp2 <- subset(temp2, Stability.index>0.8)$GeneSymbol
temp3 <- subset(temp3, Stability.index>0.8)$GeneSymbol


temp4 = intersect(temp, temp2)
temp5 = intersect(temp3,temp4)

myHK <- temp5
length(myHK)
#[1] 26


HKsearch = toupper(rowData(filtered)$Symbol)

HKlist = which(HKsearch %in% myHK)

HKIDs = rowData(filtered)$ID[HKlist]

HKsymbols = rowData(filtered)$Symbol[HKlist]





#-------------------------------------------------------------------------------------------------
#getting logcounts and deconv factors, making them sce objects

set.seed(100)
#I already have clust.filtered
filteredsum <- computeSumFactors(filtered, cluster=clust.filtered, min.mean=0.1)
libsf.sce <- logNormCounts(filteredsum, size_factors=libsf.filtered)
deconv.sce <- logNormCounts(filteredsum, size_factors=deconv.filtered)
assayNames(filteredsum)

save(deconv.sce,file="Results/Deconvolution/deconv.sce.RData",compress=TRUE)
save(libsf.sce,file="Results/LibrarySize/libsf.sce.RData",compress = TRUE)

#--------------------------------------------------------------------------------------------
# #finding the key for our housekeeping genes
# 
# rowData(filtered)$Symbol <- toupper(rowData(filtered)$Symbol)
# which("SRP14"==rowData(filtered)$Symbol)
# #3707
# which("GDI2"==rowData(filtered)$Symbol)
# #24890
# which("AHSA1"==rowData(filtered)$Symbol)
# #26780
# which("RBX1"==rowData(filtered)$Symbol)
# #28008
# which("SON"==rowData(filtered)$Symbol)
# #29379
# which("CSNK2B"==rowData(filtered)$Symbol)
# #30178
# 
# rownames(filtered)[3707]
# #ENSMUSG00000009549
# rownames(filtered)[24890]
# #ENSMUSG00000021218
# rownames(filtered)[26780]
# #ENSMUSG00000021037
# rownames(filtered)[28008]
# #ENSMUSG00000022400
# rownames(filtered)[29379]
# #ENSMUSG00000022961
# rownames(filtered)[30178]
# #ENSMUSG00000024387

#---------------------------------------------------------------------
#graphing after choosing standard

library(scater)

#testing for library size factors

#9 at a time
gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
#                 colour_by="time",
                 features="ENSMUSG00000006412") +
  ggtitle("Library Size Factors (1)"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000009549"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000074884"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000027620"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000002015"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000006699"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000029038"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000038803"),
plotExpression(dimred.sce,
               x="label", 
               #                 colour_by="time",
               features="ENSMUSG00000015804")
)

dev.copy2pdf(file="Results/HKLib1.pdf")

gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000039771") +
    ggtitle("Library Size Factors (2)"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000019494"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000051695"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000030298"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000004667"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000074781"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000031879"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000084786"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000020457")
)

dev.copy2pdf(file="Results/HKLib2.pdf")

gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000021218") +
    ggtitle("Library Size Factors (3)"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000021037"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000022427"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000022400"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000001289"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000022961"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000024387"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000040385"),
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000096994")
)

dev.copy2pdf(file="Results/HKLib3.pdf")

#dev.copy2pdf(file="Data/EB/after-normalization.pdf")

#do this for both types of normalization

#2 normalizations, do 1, tell why
#normalization by library size
#sce analysis


#need to classify samples to go ahead

#pick 1 normalization method randomly to go ahead



#--------------------------------------------------------------------------
#feature selection - end at 4000 genes

#I'm picking library size factors to continue this with

library(scran)
lib.sce <- modelGeneVar(libsf.sce)

# Visualizing the fit:
fit.sce <- metadata(lib.sce)
plot(fit.sce$mean, fit.sce$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)


lib.sce[order(lib.sce$bio, decreasing=TRUE),]

#--
libsf.var <- getTopHVGs(lib.sce, n=4000)
length(libsf.var)
#[1] 4000
str(libsf.var)




#--------------------------------------------------------------------------
#dimensionality reduction - run pca first, then 2 options: follow TSNE or UMAP
library(scater)
set.seed(100) # See below.
dimred.sce <- runPCA(libsf.sce, subset_row=libsf.var) 
reducedDimNames(dimred.sce)

dim(reducedDim(dimred.sce, "PCA"))

#--

# runTSNE() stores the t-SNE coordinates in the reducedDims
# for re-use across multiple plotReducedDim() calls.
set.seed(00101001101)
dimred.sce <- runTSNE(dimred.sce, dimred="PCA")
plotReducedDim(dimred.sce, dimred="TSNE")


set.seed(1100101001)
dimred.sce <- runUMAP(dimred.sce, dimred="PCA")
plotReducedDim(dimred.sce, dimred="UMAP")


#--------------------------------------------------------------------------
#Clustering - use TSNE or UMAP clusters cells - get maybe 10 groups
#then split cells by clusters for normalization

library(scran)
g <- buildSNNGraph(dimred.sce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

g.5 <- buildSNNGraph(dimred.sce, k=5, use.dimred = 'PCA')
clust.5 <- igraph::cluster_walktrap(g.5)$membership
table(clust.5)

g.15 <- buildSNNGraph(dimred.sce, k=15, use.dimred = 'PCA')
clust.15 <- igraph::cluster_walktrap(g.15)$membership
table(clust.15)

g.20 <- buildSNNGraph(dimred.sce, k=20, use.dimred = 'PCA')
clust.20 <- igraph::cluster_walktrap(g.20)$membership
table(clust.20)


library(scater)
dimred.5 <- dimred.sce
dimred.15 <- dimred.sce
dimred.20 <- dimred.sce
colLabels(dimred.20) <- factor(clust.20)
colLabels(dimred.15) <- factor(clust.15)
colLabels(dimred.sce) <- factor(clust)
colLabels(dimred.5) <- factor(clust.5)
plotReducedDim(dimred.sce, "UMAP", colour_by="label")

gridExtra::grid.arrange(
  plotReducedDim(dimred.5, "UMAP", colour_by="label"),
  plotReducedDim(dimred.sce, "UMAP", colour_by="label"),
  plotReducedDim(dimred.15, "UMAP", colour_by="label"),
  plotReducedDim(dimred.20, "UMAP", colour_by="label"),
  ncol=2
)

#dev.copy2pdf(file="Data/EB/clustering_tests_umap.pdf")


gridExtra::grid.arrange(
  plotReducedDim(dimred.5, "TSNE", colour_by="label"),
  plotReducedDim(dimred.sce, "TSNE", colour_by="label"),
  plotReducedDim(dimred.15, "TSNE", colour_by="label"),
  plotReducedDim(dimred.20, "TSNE", colour_by="label"),
  ncol=2
)

#dev.copy2pdf(file="Data/EB/clustering_tests_tsne.pdf")



#from these results, I choose 15 as my k value'
#10 was also looking good as a k value

#15 looks good for both tsne and umap


#-----------------------------------------------------------------------------
#after clustering, find expression across groups in normalization - then compare them
#classification (inside clustering)
#then, decide which method to use in normalization, follow other steps

#look for expression of marker genes from paper
#12.4

#plot paper's marker logcounts within my groups

#review what happened in the paper

#get a good number of clusters, label them via paper's markers


#marker genes for EB from paper:
#Etv2
#Gata2
#Tal1
#Zfp42
#Pou5f1
#Kdr

mymarkers = c("Etv2","Gata2","Tal1","Zfp42","Pou5f1","Kdr","Tbx20","Hand1","Mesp1")

markerids = which((rowData(filtered)$Symbol) %in% mymarkers)

markerids2 = rownames(filtered)[markerids]

markerlabels=rowData(filtered)$Symbol[markerids]

#plotExpression(dimred.sce, features=markerids2, 
#              x="label", colour_by="label")

gridExtra::grid.arrange(
  plotExpression(dimred.sce, features=markerids2[1], 
                 x="label", colour_by="label",
                 xlab=markerlabels[1]),
  plotExpression(dimred.sce, features=markerids2[2], 
                 x="label", colour_by="label",
                 xlab=markerlabels[2]),
  plotExpression(dimred.sce, features=markerids2[3], 
                 x="label", colour_by="label",
                 xlab=markerlabels[3]),
  plotExpression(dimred.sce, features=markerids2[4], 
                 x="label", colour_by="label",
                 xlab=markerlabels[4]),
  plotExpression(dimred.sce, features=markerids2[5], 
                 x="label", colour_by="label",
                 xlab=markerlabels[5]),
  plotExpression(dimred.sce, features=markerids2[6], 
                 x="label", colour_by="label",
                 xlab=markerlabels[6]),
  plotExpression(dimred.sce, features=markerids2[7], 
                 x="label", colour_by="label",
                 xlab=markerlabels[7]),
  plotExpression(dimred.sce, features=markerids2[8], 
                 x="label", colour_by="label",
                 xlab=markerlabels[8]),
  plotExpression(dimred.sce, features=markerids2[9], 
                 x="label", colour_by="label",
                 xlab=markerlabels[9]),
  ncol=3
)


dev.copy2pdf(file="Data/EB/expression_markers.pdf")
#compare this with graph from paper for presentation

which(rowData(filtered)$Symbol == "Sox17")

plotExpression(dimred.sce, features=rowData(filtered)$ID[6], 
               x="label", colour_by="label",
               xlab=markerlabels[1])


#cardiac = label 2
#hemangiogenic = label 3
#smooth muscle = label 5
#naive pluripotent = label 6
#endoderm = label 10
#primordial germ = label 13



#get 17 clusters
#look at pdgfra and kdr correspondence

#smallest p values are good indicators for genes in clusters

#look at groups and see how sure we are of clusters
#not certain about groups
#only look at named clusters, use table to confirm that we're right


#-------------------
#confirming that our cluster names are correct

#cardiac = label 2 - Lgals3 - looks like a combo of both 1 and 2
#hemangiogenic = label 3 - Cbfa2t3 - looks right
#smooth muscle = label 5 - Col1a1 - looks right
#naive pluripotent = label 6 (already pretty sure about this one)
#endoderm = label 10 - Sox17 - might be both 10 and 12 together, this was how the other marker looked
#primordial germ = label 13 - Lefty1,Dnd1 - looks more like label 11 from this

confirmations = c("Lgals3","Cbfa2t3","Col1a1","Sox17","Lefty1","Dnd1")

confirmation = which((rowData(filtered)$Symbol) %in% confirmations)

confirmations2 = rownames(filtered)[confirmation]

confirmationlabels=rowData(filtered)$Symbol[confirmation]

gridExtra::grid.arrange(
  plotExpression(dimred.sce, features=confirmations2[1], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[1]),
  plotExpression(dimred.sce, features=confirmations2[2], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[2]),
  plotExpression(dimred.sce, features=confirmations2[3], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[3]),
  plotExpression(dimred.sce, features=confirmations2[4], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[4]),
  plotExpression(dimred.sce, features=confirmations2[5], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[5]),
  plotExpression(dimred.sce, features=confirmations2[6], 
                 x="label", colour_by="label",
                 xlab=confirmationlabels[6]),
  ncol=2
)





#now, select normalization type
#compare normalization across groups, 

#save interested clusters into object
#load normalization, and replot results
#color the normalization cells according to cluster ids
#replot it, coloring by the clusters


#-----------------------------------------------------------------------------
#check doublets - hypothesis that doublets cluster
#doublet detection with clusters - hypothesis that doublet clusters are likely
# to be tipping point
#group of cells have more than 1 feature of the cluster

dimred.sce

library(scran)
dbl.out <- doubletCluster(dimred.sce)
dbl.out

library(scater)
chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$N, 
                                              type="lower", log=TRUE)]
chosen.doublet

markers <- findMarkers(dimred.sce, direction="up")
dbl.markers <- markers[[chosen.doublet]]

library(scater)
chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(dimred.sce, order_columns_by="label", features=chosen, 
            center=TRUE, symmetric=TRUE, zlim=c(-5, 5))

dbl6 = which(colData(dimred.sce)$label == 6)

dbl6genes = rowData(dimred.sce)$Symbol[dbl6]

#rp23 - clust 17
#lrrfip1 - clust 13
#arpc5 - 7,13,14
#tagln2 - 6,7,13,14

#idk what to do from here

#----------------------------------------------
#finding marker genes
topvals = 10
#adjust analysis with this

markers.binom <- findMarkers(dimred.sce, test="binom", direction="up")
names(markers.binom)

library(scater)

markerEB <- list()
for (i in 1:13){
interesting.binom <- markers.binom[[i]]
top.genes <- row.names(interesting.binom[interesting.binom$Top <= topvals,])
#top.genes <- head(rownames(interesting.binom), topvals)
barcode.list <- which(rowData(dimred.sce)$ID %in% top.genes)
markerEB[[i]] = rowData(dimred.sce)$Symbol[barcode.list]
}

#conversion from barcodes to gene names

plotExpression(dimred.sce, x="label", features=top.genes)

plotExpression(dimred.sce, x="label", features="ENSMUSG00000026124",xlab="gene")



#-----------
#copying table in

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

save(markerS1, file="Results/markerS1.RData", compress=TRUE)

#identify own markers for 13 subclusters
#translate from gene symbol

#now compare top markers for clusters with author's top markers


#-----------------------------------------
#analyzing binomial tests

#our old assumptions
#cardiac = label 2 - Lgals3 - looks like a combo of both 1 and 2
#hemangiogenic = label 3 - Cbfa2t3 - looks right
#smooth muscle = label 5 - Col1a1 - looks right
#naive pluripotent = label 6 (already pretty sure about this one)
#endoderm = label 10 - Sox17 - might be both 10 and 12 together, this was how the other marker looked
#primordial germ = label 13 - Lefty1,Dnd1 - looks more like label 11 from this

#--------
#new info

#cardiac - matches for 8
#hemangiogenic - 3 still correct
#smooth muscle - highest proportion is 5
#naive pluripotent - 6 still correct
#endoderm - 10 and 12 still correct
#primordial germ - 11 looks correct

intersects <- list()
for (i in 1:length(markerEB)){
intersects[[i]] = unlist(lapply(markerS1,function(X) intersect(X,markerEB[[i]])))
}


#figure out where cardiac is

#repeat analysis with t test

#-----------------------------------------
#t-test marker gene analysis
markers.ttest <- findMarkers(dimred.sce)
names(markers.ttest)

library(scater)
markerEBttest <- list()


for (i in 1:13){
  interesting.ttest <- markers.ttest[[i]]
  top.genes <- interesting.ttest[interesting.ttest$Top <= topvals,]
  barcode.list <- which(rowData(dimred.sce)$ID %in% row.names(top.genes))
  markerEBttest[[i]] = rowData(dimred.sce)$Symbol[barcode.list]
}

intersects.ttest <- list()
for (i in 1:length(markerEBttest)){
  intersects.ttest[[i]] = unlist(lapply(markerS1,function(X) intersect(X,markerEBttest[[i]])))
}

targets = c("cardiac", "hemangiogenic", "smooth muscle", "naive pluripotent", "endoderm", "primordial germ")

par(mfrow = c(2,3))

for(i in targets){
  temp = lapply(intersects.ttest, function(x) length(grep(i, names(x))) / length(x))
  plot(unlist(temp),main=i)
}

dev.copy2pdf(file=paste0("Results/ttestgroupproportions_",topvals,".pdf"))


#run proportional analysis for binomial tests

#---------------------------------------------------------------------------
#trajectory building
#branching point should be identification

dimred.sce

library(scater)
by.cluster <- aggregateAcrossCells(dimred.sce, ids=colData(dimred.sce)$label)
centroids <- reducedDim(by.cluster, "PCA")

dmat <- dist(centroids)
dmat <- as.matrix(dmat)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)

set.seed(1000)
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "TSNE")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)

plotTSNE(dimred.sce, colour_by="label") + 
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))



#---------------------------
.map2edges <- function(points, center, edge.ends, previous) {
  all.distances <- list()
  all.pseudo <- list()
  edge.len <- list()
  
  # Computing distance of each point from each edge.
  # Edges defined from 'center' to 'edge.ends'.
  for (i in rownames(edge.ends)) {
    edge.end <- edge.ends[i,]
    delta <- center - edge.end
    max.d <- sqrt(sum(delta^2))
    delta <- delta/max.d
    
    centered <- t(t(points) - center)
    proj <- as.numeric(centered %*% delta)
    proj <- pmax(0, pmin(proj, max.d))
    mapped <- outer(proj, delta)
    
    dist <- sqrt(rowSums((centered - mapped)^2))
    all.distances[[i]] <- dist
    all.pseudo[[i]] <- proj
    edge.len[[i]] <- max.d
  }
  
  all.distances <- do.call(cbind, all.distances)
  all.pseudo <- do.call(cbind, all.pseudo)
  chosen <- colnames(all.distances)[max.col(-all.distances)]
  
  # Flipping the distance of points to the previous node,
  # in order to enforce a directional pseudo-time.
  dist.previous <- 0
  if (!is.na(previous)) {
    on.previous <- chosen==previous
    dist.previous <- edge.len[[previous]]
    previous.proj <- dist.previous - all.pseudo[on.previous,previous,drop=FALSE]
    
    if (all(on.previous)) {
      return(list(dist=dist.previous, pseudo=list(previous.proj)))
    }
  }
  
  # Filling out the branches, where points are NA for a branch's
  # pseudo-time if they were assigned to another branch.
  output <- list()
  for (leftover in setdiff(rownames(edge.ends), previous)) {
    empty <- rep(NA_real_, nrow(points))
    if (!is.na(previous)) {
      empty[on.previous] <- previous.proj
    }
    current <- chosen==leftover
    empty[current] <- all.pseudo[current,leftover]
    output[[leftover]] <- empty
  }
  
  list(dist=dist.previous, pseudo=output)
}

originals <- reducedDim(dimred.sce, "PCA")
cluster <- colLabels(dimred.sce)
starting.cluster <- names(igraph::V(mst)[igraph::degree(mst)==1])[1]
collated <- list()

latest <- starting.cluster
parents <- NA_character_ 
progress <- list(rep(NA_real_, length(cluster)))
cumulative <- 0

while (length(latest)) {
  new.latest <- new.parents <- character(0)
  new.progress <- list()
  new.cumulative <- numeric(0)
  
  for (i in seq_along(latest)) {
    curnode <- latest[i]
    all.neighbors <- names(igraph::adjacent_vertices(mst, curnode, mode="all")[[1]])
    in.cluster <- cluster==curnode 
    
    mapped <- .map2edges(originals[in.cluster,,drop=FALSE], center=centroids[curnode,], 
                         edge.ends=centroids[all.neighbors,,drop=FALSE], previous=parents[i])
    edge.len <- mapped$dist
    pseudo <- mapped$pseudo
    
    collected.progress <- list()
    for (j in seq_along(pseudo)) {
      sofar <- progress[[i]] # yes, using 'i' here.
      sofar[in.cluster] <- pseudo[[j]] + cumulative[i]
      collected.progress[[j]] <- sofar
    }
    
    all.children <- setdiff(all.neighbors, parents[i])
    if (length(all.children)==0) {
      collated[[curnode]] <- collected.progress[[1]]
    } else {
      new.latest <- c(new.latest, all.children)
      new.parents <- c(new.parents, rep(curnode, length(all.children)))
      new.progress <- c(new.progress, collected.progress)
      new.cumulative <- c(new.cumulative, rep(cumulative[i] + edge.len, length(all.children)))
    }
  }
  
  latest <- new.latest
  parents <- new.parents
  progress <- new.progress
  cumulative <- new.cumulative
}
tscan.pseudo <- do.call(cbind, collated)

plotTSNE(dimred.sce, colour_by=I(rowMeans(tscan.pseudo, na.rm=TRUE)), text_by="label") +
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))



#---------------------------------------------------------------
#todo: run analysis again with deconvolution instead of libsf
#make sure i get the same results

#for trajectory: tell them which cluster is the beginning







#--------------------------------------------------------------------------
#bioTIP analysis







#--------------------------------------------------------------------------
#identify cell cycle genes








