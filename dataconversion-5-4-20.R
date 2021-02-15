#link for code
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#failed to allocate 106 GB vector
library(Matrix)
packageVersion("Matrix")
#version 1.2.18
matrix_dir = "C:/Users/AJ/RWork/"
barcode.path <- paste0(matrix_dir, "GSM3732839_EB_barcodes.tsv")
features.path <- paste0(matrix_dir, "GSM3732839_EB_genes.tsv")
#features = genes in population
#gene, transcript, feature are all the same thing in a row
#cells, population, library are all the same things in a column
matrix.path <- paste0(matrix_dir, "GSM3732839_EB_matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


barcode1.path <- paste0(matrix_dir, "GSM3732840_YS_barcodes.tsv")
genes1.path <- paste0(matrix_dir, "GSM3732840_YS_genes.tsv")
matrix1.path <- paste0(matrix_dir, "GSM3732840_YS_matrix.mtx")
mat1 <- readMM(file = matrix.path)
feature1.names = read.delim(genes1.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode1.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat1) = barcode.names$V1
rownames(mat1) = feature.names$V1

#dgtmatrix

A.data<-as.matrix(A.data)
write.csv2(A.data, "Adata.csv")
