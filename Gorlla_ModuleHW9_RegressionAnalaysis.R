# MODULE 9 HOMEWORK


## PROBLEM 1
library(ALL)
data(ALL) 
# finding the expression values for the GRO2 GRO2 oncogene and the GRO3 GRO3 oncogene.
GRO2_index <-grep("GRO2 GRO2 oncogene",golub.gnames[,2])
GRO3_index <-grep("GRO3 GRO3 oncogene",golub.gnames[,2])
# The expression values of these two genes- GR02, GR03
GRO2_expression <-golub[GRO2_index,]
GRO3_expression <-golub[GRO3_index,]

# (a) Find the correlation between the expression values of these two genes.
cor(GRO2_expression, GRO3_expression)
#[1] 0.7966283

# (b) Find the parametric 90% confident interval for the correlation with cor.test().
cor.test(GRO2_expression, GRO3_expression,  conf.level = 0.90)
# 90 percent confidence interval: 0.6702984 0.8780861
# p-value = 2.201e-09
# Since the p-value =  2.201e-09 is very small, there is strong evidence that the two sets of expression values are correlated.

# (c) Find the bootstrap 90% confident interval for the correlation.
nboot <- 2000 # We will resample 2000 times
boot.cor <- matrix(0,nrow=nboot, ncol = 1) #A vector to save the resampled statistics
data <- cbind(GRO2_expression, GRO3_expression) #Data set with x and y in two columns.
for (i in 1:nboot){
  dat.star <- data[sample(1:nrow(data),replace=TRUE), ] #Resample the pairs
  boot.cor[i,] <- cor(dat.star[,1], dat.star[,2]) #Correlation on resampled data
}
quantile(boot.cor[,1],c(0.05,0.95)) #Find quantiles for resampled statistics
# Hence the 90% CI is (0.5839700 0.8983656 )


# (d) Test the null hypothesis that correlation = 0.64 against the one-sided alternative that correlation > 0.64 at the α = 0.05 level. What is your conclusion? Explain you reasoning supported by the appropriate R outputs.
n<-length(GRO2_expression) #sample size n = number of pairs
T.obs<- 0.64 #Observed statistic is the correlation rho
n.perm=2000 # We will permute 2000 times
T.perm = rep(NA, n.perm) #A vector to save the permuted statistics
for(i in 1:n.perm) {
  GRO2_expression.perm = sample(GRO2_expression, n, replace=F) #permute data (x only)
  T.perm[i] = cor(GRO2_expression.perm, GRO3_expression) #Permuted statistic is the correlation
}
mean(abs(T.perm)>=abs(T.obs)) #p-value for 2-sided test
# [1] 0 is less than 0.05. Therefore, we reject the the null hypothesis that correlation = 0.64 against the one-sided alternative that correlation > 0.64 at the α = 0.05 level and conclude that the correlation is greater than 0.64


## PROBLEM 2
# (a) How many of the genes have correlation values less than negative 0.5? (Those genes are highly negatively correlated with Zyxin gene).
data(golub,package = "multest")
grep("Zyxin",golub.gnames[,2])
zyxin<-golub[2124,]
cor.values<-apply(golub,1,function(x) cor(zyxin,x))
sum(cor.values< -0.5) # genes have correlation values less than negative 0.5
# [1] 85

# (b) Find the gene names for the top five genes that are most negatively correlated with Zyxin gene.
golub.gnames[order(cor.values,decreasing=FALSE)[1:5],2]
#  top five genes that are most negatively correlated with Zyxin gene.
#[1] "Macmarcks"  
#[2]"Inducible protein mRNA"                                                                                         
#[3] "C-myb gene extracted from Human (c-myb) gene, complete primary cds, and five complete alternatively spliced cds" 
#[4]"Oncoprotein 18 (Op18) gene"                                                                                     
#[5] "54 kDa protein mRNA"   

# (c) Using the t-test, how many genes are negatively correlated with the Zyxin gene? Use a false discovery rate of 0.05. (Hint: use cor.test() to get the p-values then adjust for FDR. Notice that we want a one-sided test here.)
p_val<- apply(golub,1, function(x) cor.test(zyxin,x,alternative = "less")$p.value)
p.fdr<- p.adjust(p=p_val, method= "fdr")
sum(p.fdr<0.05)
# [1] 142

## PROBLEM 3
# (a) Is there a statistically significant linear relationship between the two genes’ expression? Use appropriate statistical analysis to make the conclusion. What proportion of the GRO3 GRO3 oncogene expression’s variation can be explained by the regression on GRO2 GRO2 oncogene expression?
reg.fit <- lm (GRO3_expression ~ GRO2_expression)
reg.fit
summary(reg.fit)
# Since the p-value =  2.201e-09 is very small than 0.05 there a statistically significant linear relationship between the two genes’expression

#(b) Test if the slope parameter is less than 0.5 at the α = 0.05 level.
confint(reg.fit, level = 0.9)
# 5 %       95 %
#  (Intercept)     -0.9428580 -0.7422600
#GRO2_expression  0.2817217  0.4346801
# The upper limit is lower than 0.5. hence, slope parameter is below 0.5 at the α = 0.05 level.

#(c) Find an 80% prediction interval for the GRO3 GRO3 oncogene expression when GRO2 GRO2 oncogene is not expressed (zero expression value).
predict(reg.fit, newdata=data.frame(GRO2_expression = 0), interval = "prediction", level = 0.8)
# fit       lwr        upr
# 1 -0.842559 -1.267563 -0.4175553
# 80% prediction interval for the GRO3 GRO3 oncogene expression when GRO2 GRO2 oncogene is not expressed (zero expression value) is -1.267563 -0.4175553


#(d) Check the regression model assumptions. Can we trust the statistical inferences from the regression fit?
shapiro.test(resid(reg.fit))
# Shapiro-Wilk test is conducted to test the normality of the residuals.
# The p-value is 0.07532, which is greater than 0.05, we accept the null hypothesis which indicates that the residuals are normally distributed
# Therefore, We can trust the statistical inferences.




## PROBLEM 4
# (a) Regress stack.loss on the other three variables. What is the fitted regression equation?
sldata<- as.data.frame(stackloss[,c('stack.loss','Air.Flow','Water.Temp','Acid.Conc.')])
lin.reg<-lm(stack.loss~Air.Flow+Water.Temp+Acid.Conc.,data=sldata)
lin.reg
# Regression equation is : 
# stack.loss = 39.9197+0.7156*Airflow+1.2953*Water.Temp-0.1521 * Acidic.Conc

# (b) Do all three variables have statistical significant effect on stack.loss? What proportion of variation in stack.loss is explained by the regression on the other three variables?
summary(lin.reg)
# For air-flow and water.temp, p-value is 5.8e-05 and 0.00263 which are less than 0.05. Air flow and water temperature have significance for stack.loss. 
# where as the p-value for Acid.conc is 0.34405  which is greater than 0.05 and does not have statistical significance 
# Therefore, of all three variables,  air-flow and water.temp have statistical significant effect on stack.loss


#(c) Find a 90% confidence interval and 90% prediction interval for stack.loss when Air.Flow=60, Water.Temp=20 and Acid.Conc.=90.
predict(lin.reg,newdata = data.frame(Air.Flow=60,Water.Temp=20,Acid.Conc.=90), interval = "confidence", level = 0.9) #90% confidence interval -(13.50069 16.96617)
# fit      lwr      upr
# 1 15.23343 13.50069 16.96617
predict(lin.reg,newdata = data.frame(Air.Flow=60,Water.Temp=20,Acid.Conc.=90), interval = "prediction", level = 0.9) #90% prediction interval ( 9.331184 21.13568)
#fit      lwr      upr
#1 15.23343 9.331184 21.13568









