library(Matrix)
matrix_dir = "/Users/ArjunKinstlick/UChicago/BioTip_Work/GSE130146_RAW/"

#loading EB
barcode.path <- paste0(matrix_dir, "GSM3732839_EB_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "GSM3732839_EB_genes.tsv.gz")
matrix.path <- paste0(matrix_dir, "GSM3732839_EB_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

#loading YS
barcode1.path <- paste0(matrix_dir, "GSM3732839_EB_barcodes.tsv.gz")
features1.path <- paste0(matrix_dir, "GSM3732839_EB_genes.tsv.gz")
matrix1.path <- paste0(matrix_dir, "GSM3732839_EB_matrix.mtx.gz")
mat1 <- readMM(file = matrix1.path)
feature1.names = read.delim(features1.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode1.names = read.delim(barcode1.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

#exporting both to csv format
#https://stackoverflow.com/questions/44554407/dgtmatrix-export-to-csv
mat.data<-as.matrix(mat)
#I ran out of memory here

write.csv(mat,"eb-matrix.csv")
write.csv(mat1,"ys-matrix.csv")

