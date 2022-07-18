# MODULE 10 HOMEWORK

## PROBLEM 1 

# Clustering analysis on the "CCND3 Cyclin D3" gene expression values of the Golub et al. (1999) data.
library(multtest)
data(golub)
clusdata <- data.frame(golub[1042,])
gol.fac <-factor(golub.cl,levels = 0:1, labels = c("ALL","AML"))

#(a) Conduct hierarchical clustering using single linkage and Ward linkage. Plot the cluster dendrogram for both fit. Get two clusters from each of the methods.
# Use function table() to compare the clusters with the two patient groups ALL/AML. Which linkage function seems to work better here?
# single linkage  
hc_single <- hclust(dist(clusdata,method="euclidean"), method = "single")
plot(hc_single,labels=gol.fac)
rect.hclust(hc_single,2)
clusters_single <- cutree(hc_single,k=2)
# Ward linkage
hc_ward <- hclust(dist(clusdata,method="euclidean"), method = "ward.D2")
plot(hc_ward,labels=gol.fac)
rect.hclust(hc_ward,2)
clusters_ward <- cutree(hc_ward,k=2)

# storing the AML and ALL labels as a dataframe
labels <- data.frame(gol.fac)
# clusters as data frame
clusters_single_df <- data.frame(clusters_single)
clusters_ward_df <- data.frame(clusters_ward)
# table for cluster results 
table(labels$gol.fac, clusters_single_df$clusters)
#     1  2
#ALL 27  0
#AML 10  1
table(labels$gol.fac, clusters_ward_df$clusters)
#    1  2
#ALL 21  6
#AML  0 11
clusters_single_df
# on cutting the dendrograms at 2 clusters: while using single linkage all the ALL were clustered together but 1o AML were misclustered(they were present in the same cluster as ALL)
# out of all the AML, 2 were present on different hierarchies each as a single cluster (had different heights due to diversification), which indicates that these AML patients are distantly related to rest of AML,
#That indicates that they are even more distant from AML than AMl are from ALL(which is false)
# While using ward linkage all AML were clustered together
#and only 6 ALL were misclustered and the first branching itself created two major clusters clearly separating the AML and ALL groups
# Therefore, ward linkage performs better.


# (b) Use k-means cluster analysis to get two clusters. Use table() to compare the two clusters with the two patient groups ALL/AML. 
clusters_km <- kmeans(clusdata,2)
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
table(labels$gol.fac, clusters_km$cluster)
#     1  2
#ALL 22  5
#AML  1 10

# (c) Which clustering approach (hierarchical versus k-means) produce the best matches to the two diagnose groups ALL/AML?
# hierarchical clustering and kmeans clustering with ward linkage has produced almost similar results where 
# kmeans misclustered 5 ALL and 1 AML patients
#     1  2
#ALL 22  5
#AML  1 10

# hierarchical clustering with ward linkage has misclustered 6 ALL patients
#    1  2
#ALL 21  6
#AML  0 11

# (d) Find the two cluster means from the k-means cluster analysis. Perform a bootstrap on the cluster means.
# Do the confidence intervals for the cluster means overlap? Which of these two cluster means is estimated more accurately?

initial <- clusters_km$centers 
n <- dim(clusdata)[1]
nboot <- 1000 
boot_cl <- matrix(NA,nrow=nboot,ncol = 2) 

for (i in 1:nboot){ 
  dat.star <- clusdata[sample(1:n,replace=TRUE),]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot_cl[i,] <- c(cl$centers[,1])}

apply(boot_cl,2,mean) 
# [1] 2.0353376 0.7032434
quantile(boot_cl[,1],c(0.025,0.975)) 
####2.5%    97.5% 
#  1.833972 2.19921

quantile(boot_cl[,2],c(0.025,0.975)) 
####2.5%    97.5% 
#  0.168568 1.071967 

# The two cluster means from the k-means cluster analysis are:
# cluster 1:  2.0353376
 # cluster 2: 0.7032434
# The two clusters remained stable over re sampling. 
# The bootstrap means that we got are are:
####2.5%    97.5% 
#  1.833972 2.19921
####2.5%    97.5% 
#  0.168568 1.071967
# which are almost same as the 2 means from the original clustering So there is less bias. 
# The estimations of cluster means are quite precise because the 95% bootstrap confidence intervals are fairly small.
# But the first bootstrap mean has been estimated more accurately as it is closer to the actual first mean as compared to the difference between 
# the second bootstrap mean and second original mean and the 95% bootstrap confidence intervals for the first mean are are also smaller.
#The difference between the bootstrap means and the k-means from the original data gives an estimate of the estimation bias. 
#It can be observed that the bias is small. The estimation is quite precise because the 95% bootstrap confidence intervals are fairly small.


# (e) Produce a plot of K versus SSE, for K=1, …, 30. How many clusters does this plot suggest?
K<-(1:30); sse<-rep(NA,length(K))
for (k in K) {
  sse[k]<-kmeans(clusdata, centers=k,nstart = 10)$tot.withinss
}
plot(K, sse, type='o', xaxt='n'); axis(1, at = K, las=2)
# According to the plot there should be 4 clusters as at that point the plot rapidly change and there is a sharp decrease in the slope of the curve, that is from this point  sse decreases slowly and becomes almost constant


## PROBLEM 2
# Cluster analysis on part of Golub data.
# (a) Select the oncogenes and antigens from the Golub data. (Hint: Use grep() ).
library(multtest);data(golub);
og <- grep('oncogene', golub.gnames[,2])
# Selecting antigens
ag <- grep('antigen', golub.gnames[,2])
# (b) table for kmeans
clusdata<-rbind(golub[og,], golub[ag,])
g.name<-rep(c("oncogene","antigen"), c(length(og), length(ag)))
k.means <- kmeans(clusdata, centers=2)
kmean.table <- table(g.name, k.means$cluster)
kmean.table
# g.name      1  2
#   antigen  34 41
#   oncogene 20 22
# table for k_mediod
library(cluster)
k_mediod<- pam(clusdata, k=2)
kmedoid.table <- table(g.name, k_mediod$cluster)
kmedoid.table
# g.name      1  2
#   antigen  49 26
#   oncogene 29 13
# both kmeans and kmedoids are not able to group antigens and oncogenes clusters perfectly

# (c) Use appropriate tests (from previous modules) to test the marginal independence in the two by two tables in (b). Which clustering method provides clusters related to the two gene groups? 
chisq.test(kmean.table)
# Pearson's Chi-squared test with Yates' continuity correction
# data:  kmean.table
# X-squared = 0.0019898, df = 1, p-value = 0.9644

chisq.test(kmedoid.table)
# Pearson's Chi-squared test with Yates' continuity correction
# data:  kmedoid.table
# X-squared = 0.041786, df = 1, p-value = 0.838

# the p-value for chi-square test on the table of k-means gives a  p-value = 0.9644, and for the table of k medoid the p-value = 0.838, 
# as both the p-values are greater than 0.05 none of the clustering method produced  clusters that are significantly related to the two gene groups 

# (d) Plot the cluster dendrograms for this part of golub data with single linkage and complete linkage, using Euclidean distance. 
par(mfrow=c(1,2))
plot(hclust(dist(clusdata, method="euclidian"), method="single"),labels=g.name)
plot(hclust(dist(clusdata, method="euclidian"), method="complete"),labels=g.name)

## PROBLEM 3

# Clustering analysis on NCI60 cancer cell line microarray data (Ross et al. 2000) 
# We use the data set in package ISLR from r-project (Not Bioconductor). You can use the following commands to load the data set. 
install.packages('ISLR')
library(ISLR)
ncidata<-NCI60$data
ncilabs<-NCI60$labs
# The ncidata (64 by 6830 matrix) contains 6830 gene expression measurements on 64 cancer cell lines. The cancer cell lines labels are contained in ncilabs. We do clustering analysis on the 64 cell lines (the rows). 
# (a) Using k-means clustering, produce a plot of K versus SSE, for K=1,…, 30. How many clusters appear to be there?

K2 <- c(1:30)
sse2 <- rep(NA,length(K2))
for(k in K2){
  sse2[k] <- kmeans(t(ncidata), centers = k, nstart = 10)$tot.withinss
}
plot(K2, sse2, type = 'o', xaxt = 'n')
axis(1, at = K2, las = 2)
# from the plot there appears to be 4 or 5 clusters

# (b) Do K-medoids clustering (K=7) with 1-correlation as the dissimilarity measure on the data.
# Compare the clusters with the cell lines. Which type of cancer is well identified in a cluster? 
# Which type of cancer is not grouped into a cluster? According to the clustering results, which types of cancer are most similar to ovarian cancer?
# For (b) make sure you show the table in the output file based on which you are making these conclusions.
cl.7medoids <- pam(as.dist(1-cor(t(ncidata))), k = 7)
table(ncilabs,cl.7medoids$clustering)
# ncilabs       1 2 3 4 5 6 7
#   BREAST      0 3 0 0 2 0 2
#   CNS         1 4 0 0 0 0 0
#   COLON       0 0 0 7 0 0 0
#   K562A-repro 0 0 0 0 0 1 0
#   K562B-repro 0 0 0 0 0 1 0
#   LEUKEMIA    0 0 0 0 0 6 0
#   MCF7A-repro 0 0 0 0 1 0 0
#   MCF7D-repro 0 0 0 0 1 0 0
#   MELANOMA    0 1 0 0 0 0 7
#   NSCLC       2 2 0 3 1 1 0
#   OVARIAN     2 0 1 2 1 0 0
#   PROSTATE    0 0 1 1 0 0 0
#   RENAL       7 1 1 0 0 0 0
#   UNKNOWN     0 0 1 0 0 0 0

# Melanoma is well identified in a cluster as all the 7 Melanoma cell lines , are only present in a single cluster (7th cluster) with only one dissimilar
# #   MELANOMA    0 1 0 0 0 0 7

# NSCLC is not grouped into a cluster as it's cell lines are all distributed among 5 different clusters
# #   NSCLC       2 2 0 3 1 1 0

# according to the clustering results, the types of cancer are most similar to ovarian cancer are  NSCLC,renal, colon

#  # NSCLC       2 2 0 3 1 1 0
# # RENAL        7 1 1 0 0 0 0
# #  COLON       0 0 0 7 0 0 0
