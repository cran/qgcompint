## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----betatab, include=TRUE, echo=FALSE----------------------------------------
bn = c(
  paste0("psi_", 0),
  paste0("beta_", 1:7),
  paste0("psi_", 2),
  paste0("eta_", 1:7))
bv = c(0,
  c(0.8,0.6,0.3,-0.3,-0.3,-0.3, 0),
  0,
  c(1.0,0.0,0.0,0.0,0.2,0.2,0.2))

dt = data.frame(value=bv, row.names = bn)
print(dt)

## ----intro_gen, include=TRUE--------------------------------------------------
 library(qgcompint)
 set.seed(42)
 dat1 <- simdata_quantized_emm(
  outcometype="continuous",
# sample size
  n = 300,
# correlation between x1 and x2,x3,...
  corr=c(0.8,0.6,0.3,-0.3,-0.3,-0.3),    
# model intercept
  b0=0,
# linear model coefficients for x1,x2,... at referent level of interacting variable
  mainterms=c(0.3,-0.1,0.1,0.0,0.3,0.1,0.1), 
# linear model coefficients for product terms between x1,x2,... and interacting variable  
  prodterms = c(1.0,0.0,0.0,0.0,0.2,0.2,0.2),
# type of interacting variable
  ztype = "binary",                        
# number of levels of exposure
  q = 4,                                   
# residual variance of y
  yscale = 2.0                            
)
names(dat1)[which(names(dat1)=="z")] = "M"

print("data")
head(dat1)
print("modifier")
table(dat1$M)
print("outcome")
summary(dat1$y)
print("exposure correlation")
cor(dat1[,paste0("x",1:7)])

## ----first_step_fit, include=TRUE---------------------------------------------
qfit1 <- qgcomp.emm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat1,
  expnms = paste0("x",1:7),
  emmvar = "M",
  q = 4)
qfit1

## ----first_step_bound, include=TRUE-------------------------------------------
pointwisebound(qfit1, emmval=0)
pointwisebound(qfit1, emmval=1)

## ----first_step_plot, include=TRUE, fig.width=5,fig.height=4, fig.cap=paste("Weights for M =", 0:1)----
plot(qfit1, emmval=0)
plot(qfit1, emmval=1)

## ----first_step_boot, include=TRUE--------------------------------------------
qfit1b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data=dat1,
  expnms = paste0("x",1:7),
  emmvar = "M",
  q = 4)
qfit1b

## ----first_step_boot_plot, include=TRUE, fig.width=6,fig.height=4, fig.cap=paste("Pointwise comparisons for M =", 0:1)----
plot(qfit1b, emmval=0)
plot(qfit1b, emmval=1)

## ----first_step_boot_plot_scale, include=TRUE, fig.width=6,fig.height=4, fig.cap=paste("Pointwise comparisons (same scale) for M =", 0:1)----
p1 <- plot(qfit1b, emmval=0, suppressprint = TRUE)
p2 <- plot(qfit1b, emmval=1, suppressprint = TRUE)
p1 + ggplot2::coord_cartesian(ylim=c(0,10))
p2 + ggplot2::coord_cartesian(ylim=c(0,10))

## ----catmod_gen, include=TRUE-------------------------------------------------
 set.seed(23)
 dat2 <- simdata_quantized_emm(
  outcometype="logistic",
# sample size
  n = 300,
# correlation between x1 and x2,x3,...
  corr=c(0.6,0.5,0.3,-0.3,-0.3,0.0),    
# model intercept
  b0=-2,
# linear model coefficients for x1,x2,... at referent level of interacting variable
  mainterms=c(0.1,-0.1,0.1,0.0,0.1,0.1,0.1), 
# linear model coefficients for product terms between x1,x2,... and interacting variable  
  prodterms = c(0.2,0.0,0.0,0.0,0.2,-0.2,0.2),
# type of interacting variable
  ztype = "categorical",                        
# number of levels of exposure
  q = 4,                                   
# residual variance of y
  yscale = 2.0                            
)

print("data")
head(dat2)
print("modifier")
table(dat2$z)
print("outcome")
table(dat2$y)
print("exposure correlation")
cor(dat2[,paste0("x",1:7)])

## ----cat_mod_fit_wrong, include=TRUE------------------------------------------
qfit.wrong <- qgcomp.emm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x",1:7),
  emmvar = "z",
  q = 4, family=binomial())
qfit.wrong

## ----cat_mod_fit, include=TRUE------------------------------------------------
dat2$zfactor = as.factor(dat2$z)
# using asymptotic-based confidence intervals
qfit2 <- qgcomp.emm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x",1:7),
  emmvar = "zfactor",
  q = 4, family=binomial())
# using bootstrap based confidence intervals (estimate a)
set.seed(12312)
qfit2b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x",1:7),
  emmvar = "zfactor",
  q = 4, family=binomial(), rr = FALSE)
qfit2
qfit2b

## ----cat_mod_bound, include=TRUE,fig.width=5,fig.height=4---------------------
print("output the weights at Z=0")
getstratweights(qfit2, emmval=0)
print("output pointwise comparisons at Z=0")
pointwisebound(qfit2, emmval=0)
print("plot weights at Z=0")
plot(qfit2, emmval=0)

print("output stratum specific joint effect estimate for the mixture at Z=2")
print(getstrateffects(qfit2, emmval=2))
print("output the weights at Z=2")
print(getstratweights(qfit2, emmval=2))
print("output pointwise comparisons at Z=2")
pointwisebound(qfit2, emmval=2)
plot(qfit2, emmval=2)

print("output stratum specific joint effect estimate for the mixture at Z=2 from bootstrapped fit")
print(getstrateffects(qfit2b, emmval=2))
print("output pointwise comparisons at Z=2 from bootstrapped fit")
print(pointwisebound(qfit2b, emmval=2))
print("output modelwise confidence bounds at Z=2 from bootstrapped fit")
print(modelbound(qfit2b, emmval=2))

print("Plot pointwise comparisons at Z=2 from bootstrapped fit")
plot(qfit2b, emmval=2)


## ----contmod_gen, include=TRUE------------------------------------------------
 set.seed(23)
 dat3 <- simdata_quantized_emm(
  outcometype="continuous",
# sample size
  n = 100,
# correlation between x1 and x2,x3,...
  corr=c(0.8,0.6,0.3,-0.3,-0.3,-0.3),    
# model intercept
  b0=-2,
# linear model coefficients for x1,x2,... at referent level of interacting variable
  mainterms=c(0.3,-0.1,0.1,0.0,0.3,0.1,0.1), 
# linear model coefficients for product terms between x1,x2,... and interacting variable  
  prodterms = c(1.0,0.0,0.0,0.0,0.2,0.2,0.2),
# type of interacting variable
  ztype = "continuous",                        
# number of levels of exposure
  q = 4,                                   
# residual variance of y
  yscale = 2.0                            
)
names(dat3)[which(names(dat3)=="z")] = "CoM"

head(dat3)
summary(dat3$CoM)
summary(dat3$y)
cor(dat3[,paste0("x",1:7)])

## ----cont_mod_fit, include=TRUE-----------------------------------------------
qfit3 <- qgcomp.emm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x",1:7),
  emmvar = "CoM",
  q = 4)
qfit3
qfit3b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x",1:7),
  emmvar = "CoM",
  q = 4)
qfit3b

## ----cont_mod_bound, include=TRUE,fig.width=5,fig.height=4--------------------
print("output/plot the weights at CoM=0")
getstratweights(qfit3, emmval=0)
plot(qfit3, emmval=0)

print("output stratum specific joint effect estimate for the mixture at CoM=0")
print(getstrateffects(qfit3, emmval=0))

print("output pointwise comparisons at CoM=0")
print(pointwisebound(qfit3, emmval=0))


print("output/plot the weights at the 80%ile of CoM")
getstratweights(qfit3, emmval=quantile(dat3$CoM, .8))
plot(qfit3, emmval=quantile(dat3$CoM, .8))


print("output stratum specific joint effect estimate for the mixture at the 80%ile of CoM")
print(getstrateffects(qfit3, emmval=quantile(dat3$CoM, .8)))

print("output pointwise comparisons at at the 80%ile of CoM")
print(pointwisebound(qfit3, emmval=quantile(dat3$CoM, .8)))


## ----cont_mod_bound_boot, include=TRUE,fig.width=5,fig.height=4---------------
print("plot the pointwise effects at CoM=0")
plot(qfit3b, emmval=0)

print("output stratum specific joint effect estimate for the mixture at CoM=0")
print(getstrateffects(qfit3b, emmval=0))

print("output pointwise comparisons at CoM=0")
print(pointwisebound(qfit3b, emmval=0))


print("plot the pointwise effects at the 80%ile of CoM")
plot(qfit3b, emmval=quantile(dat3$CoM, .8))


print("output stratum specific joint effect estimate for the mixture at the 80%ile of CoM")
print(getstrateffects(qfit3b, emmval=quantile(dat3$CoM, .8)))

print("output pointwise comparisons at at the 80%ile of CoM")
print(pointwisebound(qfit3b, emmval=quantile(dat3$CoM, .8)))


## ----cont_mod_nlfit, include=TRUE---------------------------------------------
qfit3bnl <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7 + x1*(x2 + x3 + x4 + x5 + x6 +x7),
  data = dat3,
  expnms = paste0("x",1:7),
  emmvar = "CoM",
  q = 8, degree= 2)
qfit3bnl

## ----cont_mod_nl_plot_base, include=TRUE, fig.width=6,fig.height=4------------
print(pointwisebound(qfit3bnl, emmval=-1))
print(pointwisebound(qfit3bnl, emmval=-1))
print(pointwisebound(qfit3bnl, emmval=1))
plot(qfit3bnl, emmval=-1, modelband=TRUE, pointwiseref=4)
plot(qfit3bnl, emmval=1, modelband=TRUE, pointwiseref=4)

## ----cont_mod_nl_plot, include=TRUE, fig.width=6,fig.height=4-----------------
library(ggplot2)
plotdata = data=data.frame(q=qfit3bnl$index, ey=qfit3bnl$y.expected, modifier=qfit3bnl$emmvar.msm)
ggplot() + 
         geom_point(aes(x=q, y=ey, color=modifier), data=plotdata) + 
         geom_point(aes(x=q, y=ey), color="purple", data=plotdata[plotdata$modifier>1,], pch=1, cex=3) + 
         geom_smooth(aes(x=q, y=ey), se=FALSE, color="purple", data=plotdata[plotdata$modifier>1,], method = 'loess', formula='y ~ x') + 
         geom_smooth(aes(x=q, y=ey), se=FALSE, color="red", data=plotdata[plotdata$modifier < -1,], method = 'loess', formula='y ~ x') + 
         geom_point(aes(x=q, y=ey), color="red", data=plotdata[plotdata$modifier < -1,], pch=1, cex=3) + 
  theme_classic() + 
  labs(y="Expected outcome", x="Quantile score value (0 to q-1)") + 
  scale_color_continuous(name="Value\nof\nmodifier")
         


## ----survmod_gen, include=TRUE------------------------------------------------
 set.seed(23)
 dat4 <- simdata_quantized_emm(
  outcometype="survival",
# sample size
  n = 200,
# correlation between x1 and x2,x3,...
  corr=c(0.8,0.6,0.3,-0.3,-0.3,-0.3),    
# model intercept
  b0=-2,
# linear model coefficients for x1,x2,... at referent level of interacting variable
  mainterms=c(0.0,-0.1,0.1,0.0,0.3,0.1,0.1), 
# linear model coefficients for product terms between x1,x2,... and interacting variable  
  prodterms = c(0.1,0.0,0.0,0.0,-0.2,-0.2,-0.2),
# type of interacting variable
  ztype = "categorical",                        
# number of levels of exposure
  q = 4,                                   
# residual variance of y
  yscale = 2.0                            
)
dat4$zfactor = as.factor(dat4$z)
head(dat4)
summary(dat4$zfactor)
summary(dat4$time)
table(dat4$d) # 30 censored

cor(dat4[,paste0("x",1:7)])

## ----survival, include=TRUE---------------------------------------------------
qfit4 <- qgcomp.emm.cox.noboot(survival::Surv(time, d)~x1+x2+x3+x4+x5+x6+x7,
  data = dat4,
  expnms = paste0("x",1:7),
  emmvar = "zfactor",
  q = 4)

qfit4
plot(qfit4, emmval=0)
getstratweights(qfit4, emmval=2)
getstrateffects(qfit4, emmval=2)
pointwisebound(qfit4, emmval=1)

