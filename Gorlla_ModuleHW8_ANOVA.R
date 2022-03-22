# MODULE 8 HOMEWORK

# PROBLEM 1
# On the ALL data set, consider the ANOVA on the gene with the probe “109_at” expression values on B-cell patients in 5 groups: B, B1, B2, B3 and B4. 
# (a) Conduct the one-way ANOVA. Do the disease stages affect the mean gene expression value? 
# (b) From the linear model fits, find the mean gene expression value among B3 patients. Make sure you show the summary table in your submission.
# (c) Which group’s mean gene expression value is different from that of group B?
# (d) Use the pairwise comparisons at FDR=0.05 to find which group means are different. Show the output of your code. What is your conclusion?
# (e) Check the ANOVA model assumptions with diagnostic tests? Do we need to apply robust ANOVA tests here? If yes, apply the appropriate tests and state your conclusion.

## PROBLEM 1(a)
#load the package and the data set.
library(ALL)
data(ALL) 
ALLB01234 <- ALL[,ALL$BT %in% c("B", "B1","B2","B3", "B4")] #patients in 5 stages
y <- exprs(ALLB01234)["109_at",] #exprs function gives gene expression values
anova(lm(y ~ ALLB01234$BT)) #anova. Group indicator in ALLB01234$BT
# From the ANOVA table, p-value= 0.01082 is  small than 0.05 and we reject the null hypothesis. 
# Hence we conclude that in gene “109_at” the disease stages B, B1, B2, B3 and B4 affect the mean gene expression value


## PROBLEM 1(b)
summary(lm(y ~ ALLB01234$BT))  # summary of the anova results Without Intercept Term
# the B3 group mean is 6.6853 from the second output, and is calculated as
# 6.8102 -0.1249 =6.6853


## PROBLEM 1(c)
# Determining Which Group Differ
# using summary 
summary(lm(y ~ ALLB01234$BT))
# using pairwise comparison 
pairwise.t.test(y,  ALLB01234$BT)
# since the p values for the differences of all the other group means with group B are greater than 0.05 
# there is no difference between the mean of Group B and all the other groups


## PROBLEM 1(d)
pairwise.t.test(y, ALLB01234$BT,p.adjust.method='fdr')
# After the comparison we see that there is only difference between the group means of disease STAGES  B2 and B4 as the  p value is 0.01 which is less than 0.05 and are statistically significant 



## PROBLEM 1(e)
# Testing homoscedasticity
install.packages("lmtest")
library(lmtest)

# Shapiro-Wilk test to check normality of data.
shapiro.test(residuals(lm(y ~ ALLB01234$BT)))
# Since the  p-value =  0.11771 is greater than 0.05, we accept the null-hypothesis of normally distributed residuals. 
#Therefore, the normality assumption holds true

# The homoscedasticity assumption can be tested by the Breusch and Pagan test
bptest(lm(y ~ ALLB01234$BT), studentize = FALSE) #test equal variances
# Since the p-value  p-value = 0.883 is greater than 0.05 , we accept the null-hypothesis of equal variances.
#Therefore, the equal variances assumption  holds true. And assumptions is not violated
# According to the test results of  Shapiro test and Breusch-Pagan test  and the anova assumptions are not violated and there is no need to apply robust anova tests




 # PROBLEM -2
# Apply the nonparametric Kruskal-Wallis tests for every gene on the B-cell ALL patients in stage B, B1, B2, B3, B4 from the ALL data. 
# (Hint: use the apply() function.)
# (a) Use FDR adjustments at 0.05 level. How many genes are expressed different in some of the groups? 
# (b) Find the probe names for the top five genes with smallest p-values. 

## PROBLEM 2 (a)
data(ALL,package="ALL")
library(ALL)
ALLB01234 <- ALL[,ALL$BT %in% c("B", "B1","B2","B3", "B4")] #patients in 5 stages
all_genes<- exprs(ALLB01234)
pvalues <- apply(all_genes,1,function(x) kruskal.test(x ~ ALLB01234$BT)$p.value)
p.adj <- p.adjust(pvalues, method = "fdr")
# How many genes are expressed different in some of the groups
sum(p.adj< 0.05)
# [1] 423

## PROBLEM 2 (b)
# Find the probe names for the top five genes with smallest p-values. 
names(sort(p.adj)[1:5])
# [1] "1389_at"   "38555_at"  "40268_at"  "1866_g_at" "40155_at" 


                 
                 
# PROBLEM 3
# On the ALL data set, we consider the ANOVA on the gene with the probe “38555_at” expression values on two factors. The first factor is the disease stages: B1, B2, B3 and B4 (we only take patients from those four stages). The second factor is the gender of the patient (stored in the variable ALL$sex). 
# (a) Conduct the appropriate ANOVA analysis. Does any of the two factors affects the gene expression values? Are there interaction between the two factors?
# (b) Check the ANOVA model assumption with diagnostic tests? Are any of the assumptions violated?
                 
## PROBLEM 3 (a)

library("ALL")
data(ALL)
ALL_stages_gender <- ALL[,which(ALL$BT %in% c( "B1","B2","B3", "B4") & ALL$sex %in% c("M","F"))] #select patients
x<-exprs(ALL_stages_gender)["38555_at",] #gene 38555_at expression values
Bcell<-ALL_stages_gender$BT # B-cell stages
gender<-ALL_stages_gender$sex # gender 
anova(lm(x~ Bcell*gender))
# We can see that the disease stage affects the gene expression values wih a p value of 1.818e-09 
# and the gender does not affect the gene expression value as it has p-value =  0.7851   which is greater than 0.05
# and there is no statistical significant interaction between the two factors(disease stage and gender) with the  p value 0.9095

## PROBLEM 3 (b)
# Shapiro-Wilk test to check normality of data.
shapiro.test(residuals(lm(x~ Bcell*gender)))
# Since the  p-value =  0.03291 is very small than 0.05, we reject the null-hypothesis of normally distributed residuals. 
#Therefore, the normality assumption does not hold and is voilated. And we need other more robust tests
# The homoscedasticity assumption can be tested by the Breusch and Pagan test
bptest(lm(x~ Bcell*gender), studentize = FALSE) #test equal variances
# Since the p-value  p-value = 0.4539 is greater than 0.05 , we accept the null-hypothesis of equal variances. Therefore, the equal variances assumption  holds true. And assumptions is not violated



