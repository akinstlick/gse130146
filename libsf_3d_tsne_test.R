#testing with 3d TSNE

setwd("D:/SingleCellAnalysis/GSE130146/Results/LibrarySize/")

#-----------------------------------
#normalization

#the author did not provide details for their normalization, so we used normalization by deconvolution
# and normalization by library size
#-----------------------------------------
#this document uses library size to continue on later analysis
#------------------------------------------

library(scater)
libsf.filtered <- librarySizeFactors(filtered)
summary(libsf.filtered)

hist(log10(libsf.filtered), breaks=100, xlab="Log10[Size factor]", col='grey80')

library(scran)
set.seed(100)
clust.filtered <- quickCluster(filtered) 
table(clust.filtered)

head(rowData(filtered))

#----------------------------------------------------
#housekeeping

load(file="../myHK.RData")

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
assayNames(filteredsum)

save(libsf.sce,file="libsf.sce.RData",compress = TRUE)


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
dimred.sce <- runTSNE(dimred.sce, dimred="PCA",ncomponents=3)
plotReducedDim(dimred.sce, dimred="TSNE",ncomponents=3)

#I use 3 dimensions in my TSNE reduction here
#to make sure results are consistent


#--------------------------------------------------------------------------
#Clustering - use TSNE or UMAP clusters cells - get maybe 10 groups
#then split cells by clusters for normalization

library(scran)
g <- buildSNNGraph(dimred.sce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

#g.5 <- buildSNNGraph(dimred.sce, k=5, use.dimred = 'PCA')
#clust.5 <- igraph::cluster_walktrap(g.5)$membership
#table(clust.5)

#g.15 <- buildSNNGraph(dimred.sce, k=15, use.dimred = 'PCA')
#clust.15 <- igraph::cluster_walktrap(g.15)$membership
#table(clust.15)

#g.20 <- buildSNNGraph(dimred.sce, k=20, use.dimred = 'PCA')
#clust.20 <- igraph::cluster_walktrap(g.20)$membership
#table(clust.20)


library(scater)
#dimred.5 <- dimred.sce
#dimred.15 <- dimred.sce
#dimred.20 <- dimred.sce
#colLabels(dimred.20) <- factor(clust.20)
#colLabels(dimred.15) <- factor(clust.15)
colLabels(dimred.sce) <- factor(clust)
#colLabels(dimred.5) <- factor(clust.5)
#plotReducedDim(dimred.sce, "UMAP", colour_by="label")

#gridExtra::grid.arrange(
#  plotReducedDim(dimred.5, "UMAP", colour_by="label"),
#  plotReducedDim(dimred.sce, "UMAP", colour_by="label"),
#  plotReducedDim(dimred.15, "UMAP", colour_by="label"),
#  plotReducedDim(dimred.20, "UMAP", colour_by="label"),
#  ncol=2
#)

#dev.copy2pdf(file="Data/EB/clustering_tests_umap.pdf")


#gridExtra::grid.arrange(
#  plotReducedDim(dimred.5, "TSNE", colour_by="label"),
plotReducedDim(dimred.sce, "TSNE", colour_by="label",ncomponents=3)#,

#6 = progenitor
#8 = cardiac
#1 and 11 = primordial germ

dev.copy2pdf(file="libsf_tsne3d.pdf")


#  plotReducedDim(dimred.15, "TSNE", colour_by="label"),
#  plotReducedDim(dimred.20, "TSNE", colour_by="label"),
#  ncol=2
#)

#dev.copy2pdf(file="Data/EB/clustering_tests_tsne.pdf")



#from these results, I choose 15 as my k value'
#10 was also looking good as a k value

#15 looks good for both tsne and umap

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

dev.copy2pdf(file="HKLib1.pdf")

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

dev.copy2pdf(file="HKLib2.pdf")

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

dev.copy2pdf(file="HKLib3.pdf")

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
#-------------------------------------------------------------
#I complete this process automatically later on,
#so this code is irrelevant
#-------------------------------------------------------------
#------------------

#confirming that our cluster names are correct

#cardiac = label 2 - Lgals3 - looks like a combo of both 1 and 2
#hemangiogenic = label 3 - Cbfa2t3 - looks right
#smooth muscle = label 5 - Col1a1 - looks right
#naive pluripotent = label 6 (already pretty sure about this one)
#endoderm = label 10 - Sox17 - might be both 10 and 12 together, this was how the other marker looked
#primordial germ = label 13 - Lefty1,Dnd1 - looks more like label 11 from this

# #confirmations = c("Lgals3","Cbfa2t3","Col1a1","Sox17","Lefty1","Dnd1")
# 
# #confirmation = which((rowData(filtered)$Symbol) %in% confirmations)
# 
# #confirmations2 = rownames(filtered)[confirmation]
# 
# #confirmationlabels=rowData(filtered)$Symbol[confirmation]
# 
# #gridExtra::grid.arrange(
# #  plotExpression(dimred.sce, features=confirmations2[1], 
# #                 x="label", colour_by="label",
# #                 xlab=confirmationlabels[1]),
# #  plotExpression(dimred.sce, features=confirmations2[2], 
# #                 x="label", colour_by="label",
# #                 xlab=confirmationlabels[2]),
# #  plotExpression(dimred.sce, features=confirmations2[3], 
# #                  x="label", colour_by="label",
# #                  xlab=confirmationlabels[3]),
# #   plotExpression(dimred.sce, features=confirmations2[4], 
# #                  x="label", colour_by="label",
# #                  xlab=confirmationlabels[4]),
# #   plotExpression(dimred.sce, features=confirmations2[5], 
# #                  x="label", colour_by="label",
# #                  xlab=confirmationlabels[5]),
# #   plotExpression(dimred.sce, features=confirmations2[6], 
# #                  x="label", colour_by="label",
# #                  xlab=confirmationlabels[6]),
# #   ncol=2
# # )





#now, select normalization type
#compare normalization across groups, 

#save interested clusters into object
#load normalization, and replot results
#color the normalization cells according to cluster ids
#replot it, coloring by the clusters

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

dev.copy2pdf(file="trajectory_libsf_1.pdf")



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

dev.copy2pdf(file="trajectory_libsf_2.pdf")



#------------------------------------------------------------------
#more marker gene analysis, this time using table s1 data
#---------------------------------------
#opening table S1 markers for analysis
load(file="../markerS1.RData")
#---------------------------------
#finding marker genes
topvals = 10
#adjust analysis with this
#ALSO ADJUST IN TABLES1 CODE

markers.binom <- findMarkers(dimred.sce, test="binom", direction="up")
names(markers.binom)

library(scater)

markerEB <- list()
for (i in 1:length(markers.binom)){
  interesting.binom <- markers.binom[[i]]
  top.genes <- row.names(interesting.binom[interesting.binom$Top <= topvals,])
  #top.genes <- head(rownames(interesting.binom), topvals)
  barcode.list <- which(rowData(dimred.sce)$ID %in% top.genes)
  markerEB[[i]] = rowData(dimred.sce)$Symbol[barcode.list]
}



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

save(intersects, file="intersects_binom.RData", compress=TRUE)

targets = c("cardiac", "hemangiogenic", "smooth muscle", "naive pluripotent", "endoderm", "primordial germ")

par(mfrow = c(2,3))

for(i in targets){
  temp = lapply(intersects, function(x) length(grep(i, names(x))) / length(x))
  plot(unlist(temp),main=i)
}

dev.copy2pdf(file=paste0("binomgroupproportions_",topvals,".pdf"))


#figure out where cardiac is

#repeat analysis with t test

#-----------------------------------------
#t-test marker gene analysis
markers.ttest <- findMarkers(dimred.sce)
names(markers.ttest)

library(scater)
markerEBttest <- list()


for (i in 1:length(markers.ttest)){
  interesting.ttest <- markers.ttest[[i]]
  top.genes <- interesting.ttest[interesting.ttest$Top <= topvals,]
  barcode.list <- which(rowData(dimred.sce)$ID %in% row.names(top.genes))
  markerEBttest[[i]] = rowData(dimred.sce)$Symbol[barcode.list]
}

intersects.ttest <- list()
for (i in 1:length(markerEBttest)){
  intersects.ttest[[i]] = unlist(lapply(markerS1,function(X) intersect(X,markerEBttest[[i]])))
}

save(intersects.ttest, file="Results/LibrarySize/intersects_ttest.RData", compress=TRUE)


targets = c("cardiac", "hemangiogenic", "smooth muscle", "naive pluripotent", "endoderm", "primordial germ")

par(mfrow = c(2,3))

for(i in targets){
  temp = lapply(intersects.ttest, function(x) length(grep(i, names(x))) / length(x))
  plot(unlist(temp),main=i)
}

dev.copy2pdf(file=paste0("Results/LibrarySize/ttestgroupproportions_",topvals,".pdf"))


#run proportional analysis for binomial tests







#--------------------------------------------------------------------------
#bioTIP analysis







#--------------------------------------------------------------------------
#identify cell cycle genes








