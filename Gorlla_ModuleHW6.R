# MODULE 6 HOMEWORK

# Problem 1
## 1(a) 
# null hypothesis : The mean “H4/j gene” gene expression value in the ALL group is equal to -0.9
# and alternate hypothesis : The mean “H4/j gene” gene expression value in the ALL group is greater than -0.9. 
# Therefore,  (H0 : μ = -0.9 and H1:  = μ > -0.9)
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4_j_all <- golub[2972,gol.fac=="ALL"]
t.test(H4_j_all, mu=-0.9 , alternative="greater")
# test used : One Sample t-test
# p-value from r code : 0.01601
# conclusion : p-value = 0.01601,which is less than 0.05, so that the null hypothesis is rejected and the mean “H4/j gene” gene expression value in the ALL group is greater than -0.9.

# 1(b)
# null hypothesis : The mean “H4/j gene” gene expression value in ALL group does not differs from the mean “H4/j gene” gene expression value in the AML group.
# and alternate hypothesis : The mean “H4/j gene” gene expression value in ALL group differs from the mean “H4/j gene” gene expression value in the AML group.
# Therefore,  (H0 : μALL = μAML and H1: μALL ≠ μAML)
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4_j_all <- golub[2972,gol.fac=="ALL"]
H4_j_aml <- golub[2972,gol.fac=="AML"]
t.test(H4_j_all, H4_j_aml )
# test used : Welch Two Sample t-test
# p-value from r code :  0.1444
# conclusion : p-value =  0.1444,which is greater than 0.05, so that the null hypothesis is accepted and The mean “H4/j gene” gene expression value in ALL group does not differs from the mean “H4/j gene” gene expression value in the AML group.


# 1(c)
# null hypothesis: In the ALL group, the mean expression value for the “H4/j gene” gene is equal to the mean expression value for the “APS Prostate specific antigen” gene
# Alternate hypothesis: In the ALL group, the mean expression value for the “H4/j gene” gene is lower than the mean expression value for the “APS Prostate specific antigen” gene. 
# Therefore,  (H0 : μH4JALL = μ APSALL and H1: μH4JALL <  μAPSALL )
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4_j_all <- golub[2972,gol.fac=="ALL"]
APS_Psa_all<- golub[2989,gol.fac=="ALL"]
t.test(H4_j_all,APS_Psa_all, mu=0 , alternative="less", paired = T)
# test used : Paired t-test
# p-value from r code :   0.03886
# conclusion : p-value =   0.03886,which is less than 0.05, so that the null hypothesis is rejected and In the ALL group, the mean expression value for the “H4/j gene” gene is lower than the mean expression value for the “APS Prostate specific antigen” gene


# 1(d) 
# null hypothesis: In the ALL group, the proportion of patients for whom the “H4/j gene” expression values is greater than -0.6. We wish to show that pH4j in the ALL group is equal to 0.5.
# Alternate hypothesis: In the ALL group, the proportion of patients for whom the “H4/j gene” expression values is greater than -0.6. We wish to show that pH4j in the ALL group is less than 0.5.
# Therefore,  (H0 : pH4j = 0.5 and H1: pH4j < 0.5 )
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4_j_all <- golub[2972,gol.fac=="ALL"]
sum_H4_j_all<- sum (H4_j_all>-0.6)
binom.test(x=sum_H4_j_all, n=27, p=0.5, alternative = "less")
# test used : Exact binomial test
# p-value from r code : 0.1239
# conclusion : p-value =  0.1239 ,which is greater than 0.05, so that the null hypothesis is accepted and In the ALL group, the proportion of patients for whom the “H4/j gene” expression values is greater than -0.6. We wish to show that pH4j in the ALL group is not less than 0.5.


# 1(e) 
# null hypothesis: The proportion pH4j in the ALL group is same as the proportion pH4j in the AML group.
# Alternate hypothesis: The proportion pH4j in the ALL group differs from the proportion pH4j in the AML group.
# Therefore,  (H0 : pH4jALL = pH4jAML and H1 : pH4jALL ≠ pH4jAML )
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
H4_j_all <- golub[2972,gol.fac=="ALL"]
sum_H4_j_all<- sum (H4_j_all>-0.6)
APS_Psa_aml<- golub[2989,gol.fac=="AML"]
sum_APS_Psa_aml<- sum (APS_Psa_aml>-0.6)
X = c(sum_H4_j_all, sum_APS_Psa_aml)
N = c(27,11)
binom.test(x=X, n =N, p=0.5, alternative = "less")
# test used : Exact binomial test
# p-value from r code : 0.5
# conclusion : p-value =  0.5 ,which is greater than 0.05, so that the null hypothesis is accepted and The proportion pH4j in the ALL group is same as the proportion pH4j in the AML group.


##  PROBLEM 2
#  Suppose that the probability to reject a biological hypothesis by the results of a certain experiment is 0.05. This experiment is repeated 2000 times.
# (a) How many rejections do you expect?
# (b) What is the probability of less than 90 rejections?


p = 0.05
n = 2000
# 2(a) Numbeer of  rejections to be expected
n *p
# [1] 100
# 2(b) Numbeer of  rejections to be expected
pbinom(89, size = 2000, prob = 0.05)
# [1] 0.1400147



## PROBLEM 3 
# For testing H0: μ=3 versus HA: μ>3, we considers a new α=0.1 level test which rejects when  falls between  and . 
   #  (a) Use a Monte Carlo simulation to estimate the Type I error rate of this test when n=20. Do 10,000 simulation runs of data sets from the . Please show the R script for the simulation, and the R outputs for running the script. Provide your numerical estimate for the Type I error rate. Is this test valid (that is, is its Type I error rate same as the nominal α=0.1 level)? 
   #  (b)  Should we use this new test in practice? Why or why not?  

# 3(a)
n= 20
#generate 10000 data sets, each 10 observations from N(4,1) and store in a row
x.sim<-matrix(rnorm(10000*n, mean=3, sd =4), ncol=20)
#function for calculating t-test statistic
tstat<-function(x) (mean(x)-3)/sd(x)*sqrt(length(x))
#Calculate t-test statistic for each data set (row)
tstat.sim<-apply(x.sim,1,tstat)   
#Calculate the rejection rate
power.sim<-mean(tstat.sim < qt(0.4,df=n-1) & tstat.sim > qt(0.3,df=n-1))
#Display rejection rate (power) with its 95% CI
power.sim+c(-1,0,1)*qnorm(0.975)*sqrt(power.sim*(1-power.sim)/10000)  
# [1] 0.09373058 0.09960000 0.10546942
# yes, this is  test valid (that is, is its Type I error rate is 0.09960000 which is approximately equal to nominal level of α=0.1)

# 3(b)
# No, we should not use this new test in practice where Type I error rate usually gives the larger values than α=0.05(most commonly used) which results in more false positives.




## PROBLEM 4
# On the Golub et al. (1999) data set, do Welch two-sample t-tests to compare every gene’s expression values in ALL group versus in AML group. 
    # (a) Use Bonferroni and FDR adjustments both at 0.05 level. How many genes are differentially expressed according to these two criteria? 
   #  (b) Find the gene names for the top three strongest differentially expressed genes (i.e., minimum p-values). Hint: the gene names are stored in golub.gnames.

## 4 (a)
data(golub, package = "multtest")
gol.fac <- factor(golub.cl, levels =0:1, labels=c("ALL", "AML"))

# t-test -- p-value on each row for mean=0
p_values <- apply(golub, 1,function(x) t.test(x ~ gol.fac) $p.value)
# applying bonferroni adjustment
p_bon <- p.adjust(p=p_values, method = "bonferroni")
# applying fdr adjustment 
p_fdr <- p.adjust(p=p_values, method ="fdr") 
# genes which are expressed differentially 
sum(p_values<0.05)
# [1] 1078
# genes which are expressed differentially after applying bonferroni adjustment
sum(p_bon<0.05)
# [1] 103
# genes which are expressed differentially after applying fdr adjustment
sum(p_fdr<0.05)
# [1] 695

## 4(b)
o <- order (p_values, decreasing=FALSE)
golub.gnames[o[1:3],2]
# [1] "Zyxin"  "FAH Fumarylacetoacetate" "APLP2 Amyloid beta (A4) precursor-like protein 2"














