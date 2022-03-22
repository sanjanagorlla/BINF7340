# MODULE 7 HOMEWORK


## PROBLEM 1
# Wilcoxon two-sample tests to find the genes whose mean expression values are higher in the ALL group than in the AML group.
data(golub, package = "multtest")
gol.fac <- factor(golub.cl, levels = 0:1, labels = c("ALL","AML"))

## 1 (a) Use FDR adjustments at the 0.05 level. How many genes are expressed higher in the ALL group?
dw <- apply(golub, 1, function(x) wilcox.test(x~gol.fac, alternative = "greater") $p.value)
p.fdr <- p.adjust(p = dw, method = "fdr") 
# no of genes are expressed higher in the ALL group
sum(p.fdr < 0.05)
# [1] 407

# 1(b) Find the gene names for the top three genes with smallest p-values. Are they the same three genes with largest difference between the means in the ALL group and the AML group?
#  the gene names for the top three genes with smallest p-values
order <- order(dw,decreasing = F)
golub.gnames[order[1:3],2]
# "TCF3 Transcription factor 3 (E2A immunoglobulin enhancer binding factors E12/E47)" "Macmarcks"  "VIL2 Villin 2 (ezrin)"   

# three genes with largest difference between the means in the ALL group and the AML group
differ <- apply(golub,1, function(x) mean(x[gol.fac == "ALL"]) - mean(x[gol.fac == "AML"]))
differ.order <- order(differ,decreasing = TRUE) 
golub.gnames[differ.order[1:3],2]
# "TCL1 gene (T cell leukemia) extracted from H.sapiens mRNA for Tcell leukemia/lymphoma 1" "MB-1 gene" "GB DEF = (lambda) DNA for immunoglobin light chain"     

## No, the gene names for the top three genes with smallest p-values. Are not same as the three genes with largest difference between the means in the ALL group and the AML group?





## PROBLEM 2
# For the Golub et al. (1999) data set, apply the Shapiro-Wilk test of normality to every gene’s expression values in the AML group. How many genes do not pass 
# the test at 0.05 level with FDR adjustment? 
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
amlsh <- apply(golub[,gol.fac=="AML"], 1, function(x) shapiro.test(x)$p.value)
sum(amlsh<0.05)
p_fdr_amlsh <- p.adjust(p=amlsh, method ="fdr") 
# How many genes do not pass the test at 0.05 level with FDR adjustment? 
sum(p_fdr_amlsh<0.05)
# 225



## PROBLEM 3
# Gene "HOXA9 Homeo box A9" can cause leukemia (Golub et al., 1999). Use appropriate Wilcoxon two-sample tests to test if, for the ALL patients, the gene 
# "HOXA9 Homeo box A9" expresses at the same level as the “CD33” gene. 
data(golub, package = "multtest")
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
x <- golub[808,gol.fac=="ALL"]
y <-  golub[1391,gol.fac=="ALL"]
wilcox.test(x=y, y=x, paired = T,alternative="two.sided" )
# p-value = 0.01242
# Since the p-value ((p-value = 0.01242) is much smaller than α = 0.05, the conclusion is  gene  "HOXA9 Homeo box A9"  expresses at the same level as the “CD33” gene in All patients
# Wilcoxon signed rank test with continuity correction



## PROBLEM 4
data("UCBAdmissions")
apply(UCBAdmissions, 3, function(x) fisher.test(x))

# therefore, for each department

# DEPARTMENT A : p - value = 1.669e-05
# Since the p-value ((p-value =  1.669189e-05 ) is much smaller than α = 0.05, the conclusion is to reject the null hypothesis and conclude that the admission decision and gender for department A is independent
# DEPARTMENT B :p-value= 0.6771
# Since the p-value ((p-value =   0.6771 ) is greater than α = 0.05, the conclusion is to reject the alternate hypothesis and conclude that the admission decision and gender for department B is not independent
# DEPARTMENT C: p-value= 0.3866
# Since the p-value ((p-value =  0.3866  ) is much greater than α = 0.05, the conclusion is to reject the alternate hypothesis and conclude that the admission decision and gender for department C is not independent
# DEPARTMENT D : p-value =  0.5995
# Since the p-value ((p-value =   0.5995 ) is much greater than α = 0.05, the conclusion is to reject the alternate hypothesis and conclude that the admission decision and gender for department D is not independent
# DEPARTMENT E : p- value = 0.3604
# Since the p-value ((p-value =  0.3604 ) is much greater than α = 0.05, the conclusion is to reject the alternate hypothesis and conclude that the admission decision and gender for department E is not independent
# DEPARTMENT F : p - value =  0.5458
# Since the p-value ((p-value = 0.5458 ) is much greater than α = 0.05, the conclusion is to reject the alternate hypothesis and conclude that the admission decision and gender for department F is not independent





## PROBLEM 5
# null hypothesis: the variance in the ALL group is not smaller than the variance in the AML group
# alternate hypothesis: the variance in the ALL group is smaller than the variance in the AML group
  
CD33_gene <- golub[808,]
CD33_gene_ALL <- golub[808,gol.fac=="ALL"]
CD33_gene_AML <- golub[808,gol.fac=="AML"]

n <- length(CD33_gene)
# observed statistic
test.obs <- sd(CD33_gene_ALL)^2 / sd(CD33_gene_AML)^2
num.perm = 2000
test.perm = rep(NA, num.perm)
for(i in 1:num.perm) {
  data.perm = sample(CD33_gene, n, replace=F)
  test.perm[i] = sd(data.perm[gol.fac=="ALL"])^2 / 
    sd(data.perm[gol.fac=="AML"])^2
}
mean(test.obs >= test.perm)
# [1] 0.041
#  Since the p-value ((p-value = 0.041) is much smaller than α = 0.05, the conclusion is to reject the null hypothesis
# and accept the alternate hypothesis and  in "CD33” gene the variance in the ALL group is smaller than the variance in the AML group



