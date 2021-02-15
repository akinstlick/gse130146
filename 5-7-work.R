#Attempting to edit the code to work directly on the downloaded data

library(DropletUtils)
fnameeb <- file.path("D:/SingleCellAnalysis/GSE130146/Data/EB")
eb.pbmc <- read10xCounts(fnameeb, col.names=TRUE)
class(eb.pbmc)
#[1] "SingleCellExperiment"
#attr(,"package")
#[1] "SingleCellExperiment"

#This worked, however, I had to rename the files to 'barcodes.tsv', 'genes.tsv', and 'matrix.mtx'
#to fit the read format for read10xcounts. I have the files for EB and YS in different folders for ease of use

fnameys <- file.path("D:/SingleCellAnalysis/GSE130146/Data/EB")
ys.pbmc <- read10xCounts(fnameys, col.names=TRUE)
class(ys.pbmc)

#I will continue performing the analysis only on the EB dataset now, to ensure the
#process works. I will separate the analysis between files later.
library(scater)
rownames(eb.pbmc) <- uniquifyFeatureNames(
  rowData(eb.pbmc)$ID, rowData(eb.pbmc)$Symbol)
#this ran with no problems

#I'm worried about this library, since it seems to operate on human data. I used the mouse data
#from another assay
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(eb.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

#updates on 5/10

#I now try and perform quality control on my single cell experiment




