library(scran)
lib.sce <- modelGeneVar(deconv.sce)

# Visualizing the fit:
fit.sce <- metadata(deconv.sce)
plot(fit.sce$mean, fit.sce$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)


lib.sce[order(lib.sce$bio, decreasing=TRUE),]

#--
deconv.var <- getTopHVGs(lib.sce, n=4000)
length(deconv.var)
#[1] 4000
str(deconv.var)

#--------------------------
#dimensionality reduction - run pca first, then 2 options: follow TSNE or UMAP
library(scater)
set.seed(100) # See below.
dimred.sce <- runPCA(deconv.sce, subset_row=deconv.var) 
reducedDimNames(dimred.sce)

dim(reducedDim(dimred.sce, "PCA"))

#--------------------------------
set.seed(00101001101)
dimred.sce <- runTSNE(dimred.sce, dimred="PCA")
plotReducedDim(dimred.sce, dimred="TSNE")

#----------------------------------

g <- buildSNNGraph(dimred.sce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

#  1   2   3   4   5   6   7   8   9  10  11 
# 97 320 164 143 271 408  79  38  88  97  26


#--------------------------------
gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000006412") +
    ggtitle("Deconvolution (1)"),
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

dev.copy2pdf(file="Results/HKDeconv1.pdf")

gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000039771") +
    ggtitle("Deconvolution (2)"),
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

dev.copy2pdf(file="Results/HKDeconv2.pdf")

gridExtra::grid.arrange(
  plotExpression(dimred.sce,
                 x="label", 
                 #                 colour_by="time",
                 features="ENSMUSG00000021218") +
    ggtitle("Deconvolution (3)"),
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

dev.copy2pdf(file="Results/HKDeconv3.pdf")


