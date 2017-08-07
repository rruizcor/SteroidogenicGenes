#################################################################################################
#########################################CLEAN SCRIPT############################################
#################################################################################################

setwd("...")

##Packages and libraries
install.packages("d3heatmap")
install.packages("gplots")
install.packages("corrgram")
install.packages("corrplot")
install.packages("mclust")
library(d3heatmap)
library(gplots)
library(corrgram)
library(corrplot)
library(mclust)

# GTEX data
mydata = read.csv("xxx/GTexDataComplete.csv", header=TRUE)
# TCGA data
tcgadata = read.csv("xxx/SeteroidogenicGenesV2.csv")
########
#Gtex with organs to compare
df3 = mydata[mydata$SiteOrgan %in% c("Adrenal Gland", "Bladder", "Blood", "Brain", "Breast", "Cervix Uteri", "Colon", "Kidney", "Lung", "Muscle", "Ovary", "Pancreas", "Prostate", "Skin", "Spleen", "Testis", "Thyroid", "Uterus"),]
df3_matrix = data.matrix(df3[,7:ncol(df3)])# convert to matrix

## first heatmat
#Variables for clustering
distance = dist(df3_matrix, method = "manhattan")
cluster = hclust(distance, method = "ward")
row_distance = dist(df3_matrix, method = "manhattan")
row_cluster = hclust(row_distance, method = "ward")
col_distance = dist(t(df3_matrix), method = "manhattan")
col_cluster = hclust(col_distance, method = "ward")

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

RowSideColors = c(    # grouping row-variables into different
  rep("gray", 3),   # categories, Measurement 1-3: green
  rep("blue", 3),    # Measurement 4-6: blue
  rep("black", 4))
#heatmap
heatmap.2(df3_matrix, 
          cellnote = df3_matrix,  # same data set for cell labels
          main = "Correlation between normal mRNA and Transduceosome genes - GTEx Data", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,12),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          #breaks=col_breaks,    # enable color transition at specified limits
          #dendrogram="row",     # only draw a row dendrogram
          Rowv = as.dendrogram(row_cluster),
          Colv = as.dendrogram(col_cluster))


##############

########## Kidney Only
##GTEX RPKM values for Normal Kidney
KidneyGTEX = mydata[mydata$SiteOrgan %in% c("Kidney"),]
KidneyGTEXmatrix = KidneyGTEX[,c(7:35)] #pull numerical columns only
KidneyGTEXmatrix = scale(KidneyGTEXmatrix, center = TRUE, scale = TRUE) #Get Z scores using scale
KidneyGTEXmatrix = data.matrix(KidneyGTEXmatrix)# matrix for plot

col.pal <- RColorBrewer::brewer.pal(9, "Reds")
NormalPlot = d3heatmap(KidneyGTEXmatrix, colors = col.pal, scale = "row", dendrogram = "row", k_row = 3)


##TCGA RSEM values for Kidney Tumors
KidneyTumors = tcgadata[tcgadata$Cancer %in% c("ccRCC", "chrRCC", "papRCC"),]
KidneyTumors = KidneyTumors[complete.cases(KidneyTumors),]#to remove NA's and incomplete cases 
KidneyTumorsmatrix = KidneyTumors[,c(4:32)] #pull numerical columns only
KidneyTumorsmatrix = scale(KidneyTumorsmatrix, center = TRUE, scale = TRUE)##Get Z scores using scale
KidneyTumorsmatrix = data.matrix(KidneyTumorsmatrix) # matrix for plot

TumorsPlot = d3heatmap(KidneyTumorsmatrix, colors = col.pal, scale = "row", dendrogram = "row", k_row = 3)

NormalPlot
TumorsPlot

DataframeWithEquivalentNormalTissue = mydata[mydata$SiteOrgan %in% c("Adrenal Gland", "Bladder", "Blood", "Brain", "Breast", "Cervix Uteri", "Colon", "Kidney", "Lung", "Muscle", "Ovary", "Pancreas", "Prostate", "Skin", "Spleen", "Testis", "Thyroid", "Uterus"),]
DataframeWithEquivalentTumors = tcgadata[tcgadata$Cancer %in% c("ACC", "UC", "AML", "GBM", "Glioma", "Breast", "SqCCx", "Colon", "ccRCC", "chrRCC", "papRCC", "LungADC", "Skin-Melanoma", "SqCLung", "Sarcoma", "Serous", "Pancreas", "Cholangio", "Pheo", "Prostate", "DLBCL", "GCT-Testis", "Thyroid", "Endometrial", "Uterine-Carcinosarc"),]

DataframeWithEquivalentNormalTissue = mydata[mydata$SiteOrgan %in% c("Adrenal Gland", "Bladder", "Blood", "Brain", "Breast", "Cervix Uteri", "Colon", "Kidney", "Lung", "Muscle", "Ovary", "Pancreas", "Prostate", "Skin", "Spleen", "Testis", "Thyroid", "Uterus"),]
gtexmatrix = DataframeWithEquivalentNormalTissue[,c(7:35)]
gtexmatrix = scale(gtexmatrix, center = TRUE, scale = TRUE)
d3heatmap(gtexmatrix, colors = col.pal, scale = "row", dendrogram = "row", k_row = 3)

DataframeWithEquivalentTumors = tcgadata[tcgadata$Cancer %in% c("ACC", "UC", "AML", "GBM", "Glioma", "Breast", "SqCCx", "Colon", "ccRCC", "chrRCC", "papRCC", "LungADC", "Skin-Melanoma", "SqCLung", "Sarcoma", "Serous", "Pancreas", "Cholangio", "Pheo", "Prostate", "DLBCL", "GCT-Testis", "Thyroid", "Endometrial", "Uterine-Carcinosarc"),]
TCGATumorsmatrix = DataframeWithEquivalentTumors[,c(4:32)] #pull numerical columns only
TCGATumorsmatrix = scale(TCGATumorsmatrix, center = TRUE, scale = TRUE)
d3heatmap(TCGATumorsmatrix, colors = col.pal, scale = "row", dendrogram = "row", k_row = 3)

###I concatenated the files below to end up with the FINAL files with Z scores to compare
write.csv(gtexmatrix, "H:/Documents/Gtex data/GTEX Z SCORES.csv")
write.csv(TCGATumorsmatrix, "H:/Documents/Gtex data/TCGA Z SCORES.csv")
write.csv(DataframeWithEquivalentNormalTissue, "H:/Documents/Gtex data/GTEX DATA FRAME.csv")
write.csv(DataframeWithEquivalentTumors, "H:/Documents/Gtex data/TCGA DATA FRAME.csv")

#################################################################################################
#########################################Normalized DATA#########################################
#################################################################################################
#################################################################################################
# GTEX data
gtex = read.csv("xxx/GTEX DATA FRAME WITH Z SCORES_FINAL.csv", header = TRUE)
# TCGA data
tcga = read.csv("xxx/TCGA DATA FRAME WITH Z SCORES_FINAL.csv", header = TRUE)
#REMOVE na's
tcga = tcga[complete.cases(tcga),]#to remove NA's

#Matrix
gtexZscores = gtex[,c(7:35)]
gtexmatrix = data.matrix(gtexZscores)
tcgaZscores = tcga[,c(4:32)] #pull numerical columns only
tcgamatrix = data.matrix(tcgaZscores)

gtex2 = cor(gtexmatrix)
corrgram(gtexmatrix, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlation of Z-Scores of mRNA expression in Normal Tissues")

tcga2 = cor(tcgamatrix)
corrgram(tcgamatrix, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlation of Z-Scores of mRNA expression in TCGA Tumors")

head(round(tcga2,2))

#FUNCTION TO CALCULATE P-VALUES
# mat : is a matrix of data 
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
pvalgtex = cor.mtest(gtexmatrix)
head(pvalgtex[, 1:5])
pvaltcga = cor.mtest(tcgamatrix)
head(pvaltcga[, 1:5])


GTEXPLOT = corrplot(gtex2, method="circle", type="upper", order="hclust", tl.col="black", tl.srt=45,
                    p.mat = pvalgtex, sig.level = 0.05, addCoef.col = "black", hclust.method="centroid",
                    insig = "blank", diag=FALSE, main="Normal Tissues")

TCGAPLOT = corrplot(tcga2, method="circle", type="upper", order="hclust", tl.col="black", tl.srt=45,
                    p.mat = pvaltcga, sig.level = 0.05, addCoef.col = "black", hclust.method="centroid",
                    insig = "blank", diag=FALSE, main="Tumors")

#
col = colorRampPalette(c("blue", "white", "red"))(20)
#heatmap(x = gtexmatrix, col = col, symm = TRUE)
hc_dist= dist(gtexmatrix) 
hr_dist= dist(t(gtexmatrix))
hc_clust = hclust(hc_dist)
hr_clust = hclust(hr_dist)
heatmap.2(gtexmatrix, col=col, labRow=NA, density.info="none", scale="row",trace="none",
          Colv=as.dendrogram(hr_clust),
          Rowv=as.dendrogram(hc_clust))


#heatmap(x = tcgamatrix, col = col, symm = TRUE)
hc_dist2= dist(tcgamatrix) 
hr_dist2= dist(t(tcgamatrix))
hc_clust2 = hclust(hc_dist2)
hr_clust2 = hclust(hr_dist2)
heatmap.2(tcgamatrix, col=col, labRow=NA, density.info="none", scale="row",trace="none",
          Colv=as.dendrogram(hr_clust2),
          Rowv=as.dendrogram(hc_clust2))

#2*pnorm(-abs(gtexZscores$HSD3B2))
#2*pnorm(-abs(tcgaZscores$HSD3B2))

#plot(gtexZscores$HSD3B2)
#plot(tcgaZscores$HSD3B2)
