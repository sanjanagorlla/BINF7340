# MODULE 2 HOMEWORK


# library loading statement for "golub" data
BiocManager::install("multtest")
library(multtest)
data(golub)

#PROBLEM 1

## (a) Compute the mean expression values for every gene among “ALL” patients.
gol.fac <- factor(golub.cl, levels=0:1, labels = c("ALL", "AML"))  # factor indicating the tumor class of the patients
meanALL <- apply(golub[,gol.fac=="ALL"], 1, mean)
meanALL


## (b) Compute the mean expression values for every gene among “AML” patients.
meanAML <- apply(golub[,gol.fac=="AML"], 1, mean)
meanAML


## (c) Give the biological names of the three genes with the largest mean expression value among “ALL” patients.
order_meanAll <- order(meanALL,decreasing=TRUE)
golub.gnames[order_meanAll[1:3],2]
# output :  "GB DEF = Chromosome 1q subtelomeric sequence D1S553", "37 kD laminin receptor precursor/p40 ribosome associated protein gene", RPS14 gene (ribosomal protein S14) extracted from Human ribosomal protein S14 gene"


## (d) Give the biological names of the three genes with the largest mean expression value among “AML” patients.
order_meanAML <- order(meanAML,decreasing=TRUE)
golub.gnames[order_meanAML[1:3],2]
# output: "GB DEF = mRNA fragment for elongation factor TU (N-terminus)", "GB DEF = HLA-B null allele mRNA", "Globin, Beta"                                                
####################################################
#PROBLEM 2

data(golub)

## (a) Save the expression values of the first five genes (in the first five rows) for the AML patients in a csv file “AML5.csv”.
first_five_AML<-golub[1:5, gol.fac=="AML"]
write.csv ( first_five_AML, file = "AML5.csv" )
read.csv ( "AML5.csv" )

## (b) Save the expression values of the first five genes for the ALL patients in a plain text file “ALL5.txt”.
first_five_ALL<-golub[1:5, gol.fac=="ALL"]
write.table ( first_five_ALL, file = "ALL5.txt" )
read.table( "ALL5.txt" )

## (c) Compute the standard deviation of the expression values on the first patient, of the 100 th to 200 th genes (total 101 genes).
sd(golub[100:200, 1])
# output = 0.9174976

## (d) Compute the standard deviation of the expression values of every gene, across all patients. Find the number of genes with standard deviation greater than 1.
sd_gene<-apply(golub,1,sd)
sd_gene_1<- golub[sd_gene>1,]
dim(sd_gene_1)
# output - 123 38

## (e) Do a scatter plot of the 101th gene expressions against the 102th gene expressions, label the x-axis and the y-axis with the genes’ biological names using xlab= and ylab= control options.
x <- golub[101,]
y <-golub[102,]
plot(x,y,xlab=golub.gnames[101,2], ylab = golub.gnames[102,2] )

##################################################################################

# PROBLEM 3
# library loading statement for "ALL" data
BiocManager::install("ALL")
library(ALL)
data(ALL)

# (a) Use exprs(ALL[,ALL$BT=="B1"] to extract the gene expressions from the patients in disease stage B1. Produce one histogram of these gene expressions in the this matrix
B1<-exprs(ALL[,ALL$BT=="B1"])
hist(B1)

# (b) Compute the mean gene expressions for every gene over these B1 patients.
mean_B1 <- apply(exprs(ALL[,ALL$BT=="B1"]),1, mean)
mean_B1

# (c) Give the gene identifiers of the three genes with the largest mean.
order_B1 <- order(mean_B1,decreasing=TRUE)
mean_B1[order_B1[1:3],]
#output
# AFFX-hum_alu_at        31962_at      31957_r_at  #gene identifiers
#       13.41648        13.16671        13.15995  #means



###################################################################################
#PROBLEM 4

# (a) Find the type of the trees data object.
# library loading statement
data(trees)        # working with "trees" data
typeof(trees)      # type of trees data object
# output = "list"
View(trees)        # to view the complete data

# (b) Produce a figure with two overlaid scatter plots: Height versus Girth, Volume versus Girth(The Girth is on the x-axis)
# Do the Height plot with blue “+” symbols, and do the Volume plot with red “o” symbols. 
# You need to learn to set the ylim= control option so that all points from the two plots can all show up on the merged figure.
# Hint: you should use plot() then points() to create the overlaid two scatter plots.

a<-trees$Height          # Height
b<-trees$Girth           # Girth
c<-trees$Volume          # volume

# using  plot() to create Height versus Girth
plot(a~b,  ylim = c(0,90), pch = 3, col="blue", xlab="Girth", ylab="Height and volume", main="Girth V/S Height & Volume")

# points() to create the overlaid scatter plots: Volume versus Girth(The Girth is on the x-axis)
points(b,c,  ylim = c(0,90), pch = 1, col="red") # y lim is set as 10:90, control option so that all points from the two plots can all show up on the merged figure

# adding legend 
legend("bottomright",c("Height","volume"),col=c("blue","red"),lty=c(1,1))
