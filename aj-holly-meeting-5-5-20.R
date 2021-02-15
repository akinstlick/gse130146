#holly, aj meeting 5/5/20
#trying to convert dgtmatrix to r object SingleCellExperiment
#not enough memory to convert other matrix to csv, so trying this method
#https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html


#setwd("D:rwork")
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)
class(sce.pbmc)
#class:[1] "SingleCellExperiment"
#attr(,"package")
#[1] "SingleCellExperiment"

dim(sce.pbmc)
#dimension:[1]  33694 737280

#LINE 13 operates directly on our 
#get singlecellexperiment object, follow book, get classification, build subtype
#just follow book - quality control, normalization, feature selection, classification, trajectory building
#filtering cells might come after normalization
#try and make it readable as a matrix
#if it's readable, then save it
#trajectory building - try stream, monocole, and sc
#look for highly variable features
#pca analysis, tsne analysis, umap analysis
#dimension reduction chapter in book
#then trajectory analysis,

