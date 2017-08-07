

### adrenocortical neoplasms Affimetrix dataset

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", "plyr", "reshape2", "RColorBrewer", "scales", "grid", "factoextra", "cluster")
ipak(packages)

###

library(factoextra)
library(gplots)
library(corrplot)

####

df = read.csv("AffimetryxACA.csv")
names(df)
df$vdacavg = ((df$vdac1.1 + df$vdac1.2 + df$vdac1.3)/3)
df$slc25a4avg = ((df$slc25a4.1 + df$slc25a4.2 + df$slc25a4.3)/3)
df$dbiavg = ((df$dbi.1 + df$dbi.2 + df$dbi.3)/3)
df$acbdavg = ((df$acbd3.1 + df$acbd3.2)/2)
df$scarbavg = ((df$scarb1.1 + df$scarb1.2 + df$scarb1.3 + df$scarb1.4 + df$scarb1.5 + df$scarb1.6 + df$scarb1.7)/7)
df$fdxavg = ((df$fdx1.1 + df$fdx1.2 + df$fdx1.3)/3)
#cluster analysis
dfNum = df[,c(26, 29:31, 39, 40, 44:51)]#to creade Data frame with only numeric values for res.dist below
plot(dfNum)
my_data <- df[,c(3, 26, 29:31, 39, 40, 44:51)]#to select tissue type (by diagnosis) and genes 
my_data <- df[,c(7, 26, 29:31, 39, 40, 44:51)]#to select clinical/hormones data and genes 
rownames(my_data) = make.names(my_data[,1], unique = TRUE)#to convert first column into label column 
my_data = my_data[,c(2:15)]#To eliminate tissue type as first column and end up with matrix
# Remove any missing value (i.e, NA values for not available)
my_data <- na.omit(my_data)
# Scale variables
my_data <- scale(my_data)#scale calculates the z-scores
# View the firt 3 rows
head(my_data, n = 3)
#The classification of observations into groups, requires some methods for measuring the distance or the (dis)similarity between the observations.
res.dist <- get_dist(dfNum, stand = FALSE, method = "pearson")
fviz_dist(res.dist, order = TRUE, show_labels = TRUE, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))# PLOT
#Determine the optimal number of clusters
fviz_nbclust(my_data, kmeans, method = "gap_stat")
#Compute and visualize k-means clustering
km.res <- kmeans(my_data, 3, nstart = 25)
# Visualize
fviz_cluster(km.res, data = my_data, frame.type = "convex")+
  theme_minimal()
#PAM clustering: Partitioning Around Medoids. Robust alternative to k-means clustering, less sensitive to outliers.
# Compute PAM
pam.res <- pam(my_data, 3)
# Visualize
fviz_cluster(pam.res)
#Hierarchical clustering is an alternative approach to k-means clustering for identifying groups in the dataset. It does not require to pre-specify the number of clusters to be generated.
# 2. Compute dissimilarity matrix
d <- dist(my_data, method = "euclidean")
# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "ward.D2" )
# Enhanced k-means clustering
res.km = eclust(my_data, "kmeans", nboot = 10)
#Plot
fviz_cluster(res.km)
# Optimal number of clusters using gap statistics
res.km$nbclust
# Print result
res.km


#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################
my_data = data.matrix(my_data)# transform column 1-26 into a matrix

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

#########################################################
### D) Distance Matrix
#########################################################

distance = dist(my_data, method = "manhattan")
cluster = hclust(distance, method = "ward")
row_distance = dist(my_data, method = "manhattan")
row_cluster = hclust(row_distance, method = "ward")
col_distance = dist(t(my_data), method = "manhattan")
col_cluster = hclust(col_distance, method = "ward")

RowSideColors = c(    # grouping row-variables into different
  rep("gray", 3),   # categories, Measurement 1-3: green
  rep("blue", 3),    # Measurement 4-6: blue
  rep("black", 4))

gc()
heatmap.2(my_data,
          cellnote = my_data,  # same data set for cell labels
          #main = "Heatmap of Steroidogenic genes mRNA expression by Hormones Affimetrix Data", # heat map title
          main = "Heatmap of Steroidogenic genes mRNA expression by Diagnosis Affimetrix Data", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,12),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          #breaks=col_breaks,    # enable color transition at specified limits
          #dendrogram="row",     # only draw a row dendrogram
          Rowv = as.dendrogram(row_cluster),
          Colv = as.dendrogram(col_cluster))


##########################################################################################################
cor(my_data)

M = cor(my_data, method = "pearson")

corrplot(M, method = "circle")

#calculate means with NAs removed
plot(cyp11a1~star, data=my_data)

#use lm to fit regression line - linear regression
df = read.csv("AffimetryxACA.csv")
df$vdacavg = ((df$vdac1.1 + df$vdac1.2 + df$vdac1.3)/3)
df$slc25a4avg = ((df$slc25a4.1 + df$slc25a4.2 + df$slc25a4.3)/3)
df$dbiavg = ((df$dbi.1 + df$dbi.2 + df$dbi.3)/3)
df$acbdavg = ((df$acbd3.1 + df$acbd3.2)/2)
df$scarbavg = ((df$scarb1.1 + df$scarb1.2 + df$scarb1.3 + df$scarb1.4 + df$scarb1.5 + df$scarb1.6 + df$scarb1.7)/7)
df$fdxavg = ((df$fdx1.1 + df$fdx1.2 + df$fdx1.3)/3)
my_data <- df[,c(3, 26, 29:31, 39, 40, 44:51)]#to select tissue type (by diagnosis) and genes 

model1 = lm(cyp11a1~star, data = my_data)
model1
abline(model1, col="red")
plot(model1)
termplot(model1)
summary(model1)

