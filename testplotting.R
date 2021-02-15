temp = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Mouse_atlas.txt", sep="\t",head=T)
temp2 = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Mouse_development.txt", sep="\t",head=T)
temp3 = read.table(file="D:/SingleCellAnalysis/GSE130146/Data/Human_development.txt", sep="\t",head=T)

temp <- subset(temp, Stability.index>0.8)$GeneSymbol
temp2 <- subset(temp2, Stability.index>0.9)$GeneSymbol
temp3 <- subset(temp3, Stability.index>0.9)$GeneSymbol


temp4 = intersect(temp, temp2)
temp5 = intersect(temp3,temp)

venn(list(temp,temp2, temp3))
