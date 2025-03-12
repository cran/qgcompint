## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----betatab, include=TRUE, echo=FALSE----------------------------------------
bn = c(
  paste0("psi_", 0),
  paste0("beta_", 1:7),
  paste0("psi_", 2),
  paste0("eta_", 1:7))
bv = c(0,
  c(0.8, 0.6, 0.3, -0.3, -0.3, -0.3, 0),
  0,
  c(1.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2))

dt = data.frame(value=bv, row.names = bn)
print(dt)

## ----intro_gen, include=TRUE--------------------------------------------------
 library(qgcompint)
 set.seed(42)
 dat1 <- simdata_quantized_emm(
  outcometype="continuous",
# sample size
  n = 300,
# correlation between x1 and x2, x3, ...
  corr=c(0.8, 0.6, 0.3, -0.3, -0.3, -0.3),
# model intercept
  b0=0,
# linear model coefficients for x1, x2, ... at referent level of interacting variable
  mainterms=c(0.3, -0.1, 0.1, 0.0, 0.3, 0.1, 0.1),
# linear model coefficients for product terms between x1, x2, ... and interacting variable  
  prodterms = c(1.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2),
# type of interacting variable
  ztype = "binary",
# number of levels of exposure
  q = 4,
# residual variance of y
  yscale = 2.0                            
)
names(dat1)[which(names(dat1)=="z")] = "M"

## data
head(dat1)
## modifier
table(dat1$M)
## outcome
summary(dat1$y)
## exposure correlation
cor(dat1[, paste0("x", 1:7)])

## ----first_step_fit, include=TRUE---------------------------------------------
qfit1 <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat1,
  expnms = paste0("x", 1:7),
  emmvar = "M",
  q = 4)
qfit1

## ----first_step_bound, include=TRUE-------------------------------------------
pointwisebound(qfit1, emmval=0)
pointwisebound(qfit1, emmval=1)

## ----first_step_plot, include=TRUE, fig.width=5, fig.height=4, fig.cap=paste("Weights for M =", 0:1)----
plot(qfit1, emmval=0)
plot(qfit1, emmval=1)

## ----first_step_boot, include=TRUE--------------------------------------------
qfit1b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data=dat1,
  expnms = paste0("x", 1:7),
  emmvar = "M",
  q = 4)
qfit1b

## ----first_step_boot_plot, include=TRUE, fig.width=6, fig.height=4, fig.cap=paste("Pointwise comparisons for M =", 0:1)----
plot(qfit1b, emmval=0)
plot(qfit1b, emmval=1)

## ----first_step_boot_plot_scale, include=TRUE, fig.width=6, fig.height=4, fig.cap=paste("Pointwise comparisons (same scale) for M =", 0:1)----
library(ggplot2) # need to explicitly call ggplot
pp0b = plot(qfit1b, emmval=0, suppressprint=TRUE, geom_only=TRUE, 
             pointwisebars=TRUE, modelfitline=FALSE, modelband=FALSE)
pp1b = plot(qfit1b, emmval=1, suppressprint=TRUE, geom_only=TRUE, 
             pointwisebars=TRUE, modelfitline=FALSE, modelband=FALSE)
# relabel the plots by directly modifying the ggplot2 geometries.
# "col" is a column that sets the color
# "point" and "errorbar" are names given to ggplot2::geom_point and 
# ggplot2::geom_error bar geometries
pp0b$point$data$col = "M=0"
pp1b$point$data$col = "M=1"
pp0b$errorbar$data$col = "M=0"
pp1b$errorbar$data$col = "M=1"
# jitter along the x axis so the points do not overlap
pp0b$point$data$x = pp0b$point$data$x - 0.01 
pp1b$point$data$x = pp1b$point$data$x + 0.01 
pp0b$errorbar$data$x = pp0b$errorbar$data$x - 0.01 
pp1b$errorbar$data$x = pp1b$errorbar$data$x + 0.01 



suppressMessages(print(
  ggplot() + pp0b + pp1b + scale_colour_viridis_d(name="") + 
    labs(title="Bootstrap approach")
))



## ----first_step_ee, include=TRUE----------------------------------------------
qfit1ee <- qgcomp.emm.glm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data=dat1,
  expnms = paste0("x", 1:7),
  emmvar = "M",
  q = 4)
qfit1ee

## ----first_step_ee_plot, include=TRUE, fig.width=6, fig.height=4, fig.cap=paste("Pointwise comparisons for M =", 0:1)----
pp0ee = plot(qfit1b, emmval=0, suppressprint=TRUE, geom_only=TRUE, 
             pointwisebars=TRUE, modelfitline=FALSE, modelband=FALSE)
pp1ee = plot(qfit1b, emmval=1, suppressprint=TRUE, geom_only=TRUE, 
             pointwisebars=TRUE, modelfitline=FALSE, modelband=FALSE)
# relabel the plots by directly modifying the ggplot2 geometries.
# "col" is a column that sets the color
# "point" and "errorbar" are names given to ggplot2::geom_point and 
# ggplot2::geom_error bar geometries
pp0ee$point$data$col = "M=0"
pp1ee$point$data$col = "M=1"
pp0ee$errorbar$data$col = "M=0"
pp1ee$errorbar$data$col = "M=1"
# jitter along the x axis so the points do not overlap
pp0ee$point$data$x = pp0ee$point$data$x - 0.01 
pp1ee$point$data$x = pp1ee$point$data$x + 0.01 
pp0ee$errorbar$data$x = pp0ee$errorbar$data$x - 0.01 
pp1ee$errorbar$data$x = pp1ee$errorbar$data$x + 0.01 



suppressMessages(print(
  ggplot() + pp0ee + pp1ee + scale_colour_viridis_d(name="") + 
    labs(title="Estimating equations approach")
))

## ----catmod_gen, include=TRUE-------------------------------------------------
 set.seed(23)
 dat2 <- simdata_quantized_emm(
  outcometype="logistic",
# sample size
  n = 300,
# correlation between x1 and x2, x3, ...
  corr=c(0.6, 0.5, 0.3, -0.3, -0.3, 0.0),
# model intercept
  b0=-2,
# linear model coefficients for x1, x2, ... at referent level of interacting variable
  mainterms=c(0.1, -0.1, 0.1, 0.0, 0.1, 0.1, 0.1),
# linear model coefficients for product terms between x1, x2, ... and interacting variable  
  prodterms = c(0.2, 0.0, 0.0, 0.0, 0.2, -0.2, 0.2),
# type of interacting variable
  ztype = "categorical",
# number of levels of exposure
  q = 4,
# residual variance of y
  yscale = 2.0                            
)

## data
head(dat2)
## modifier
table(dat2$z)
## outcome
table(dat2$y)
## exposure correlation
cor(dat2[, paste0("x", 1:7)])

## ----cat_mod_fit_wrong, include=TRUE------------------------------------------
qfit.wrong <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "z",
  q = 4, family=binomial())
qfit.wrong

## ----cat_mod_fit, include=TRUE------------------------------------------------
dat2$zfactor = as.factor(dat2$z)
# using asymptotic-based confidence intervals
qfit2 <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  q = 4, family=binomial())
# using bootstrap based confidence intervals (estimate a)
set.seed(12312)
qfit2b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  B = 10, # set higher when using in real data
  q = 4, family=binomial(), rr = FALSE)
# using estimating equation based confidence intervals (estimate ee)
set.seed(12312)
qfit2ee <- qgcomp.emm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  q = 4, family=binomial(), rr = FALSE)
qfit2
qfit2b
qfit2ee

## ----cat_mod_bound, include=TRUE, fig.width=5, fig.height=4-------------------
## output the weights at Z=0
getstratweights(qfit2, emmval=0)
## output pointwise comparisons at Z=0
pointwisebound(qfit2, emmval=0)
## plot weights at Z=0
plot(qfit2, emmval=0)

## output stratum specific joint effect estimate for the mixture at Z=2
print(getstrateffects(qfit2, emmval=2))
## output the weights at Z=2
print(getstratweights(qfit2, emmval=2))
## output pointwise comparisons at Z=2
pointwisebound(qfit2, emmval=2)
plot(qfit2, emmval=2)

## output stratum specific joint effect estimate for the mixture at Z=2 from bootstrapped fit
print(getstrateffects(qfit2b, emmval=2))
## output pointwise comparisons at Z=2 from bootstrapped fit
print(pointwisebound(qfit2b, emmval=2))
## output modelwise confidence bounds at Z=2 from bootstrapped fit
print(modelbound(qfit2b, emmval=2))

## Plot pointwise comparisons at Z=2 from bootstrapped fit
plot(qfit2b, emmval=2)


## ----first_step_ee_plot_scalecat, include=TRUE, fig.width=6, fig.height=4, fig.cap=paste("Pointwise comparisons (same scale) for M =", 0:1)----
library(ggplot2) # need to explicitly call ggplot
pp0ee = plot(qfit2ee, emmval=0, suppressprint=TRUE, geom_only=TRUE, 
             modelband=TRUE, pointwisebars=FALSE, modelfitline=TRUE)
pp1ee = plot(qfit2ee, emmval=1, suppressprint=TRUE, geom_only=TRUE, 
             modelband=TRUE, pointwisebars=FALSE, modelfitline=TRUE)
pp2ee = plot(qfit2ee, emmval=2, suppressprint=TRUE, geom_only=TRUE, 
             modelband=TRUE, pointwisebars=FALSE, modelfitline=TRUE)
# relabel the plots uisn
pp0ee$line$data$col = "Z=0"
pp1ee$line$data$col = "Z=1"
pp2ee$line$data$col = "Z=2"
pp0ee$ribbon$data$col = "Z=0"
pp1ee$ribbon$data$col = "Z=1"
pp2ee$ribbon$data$col = "Z=2"
# jitter along the x axis so the points do not overlap
pp0ee$line$data$x = pp0ee$line$data$x - 0.01 
pp2ee$line$data$x = pp2ee$line$data$x + 0.01 
pp0ee$ribbon$data$x = pp0ee$ribbon$data$x - 0.01 
pp2ee$ribbon$data$x = pp2ee$ribbon$data$x + 0.01 

ggplot() + pp0ee + pp1ee + pp2ee + scale_color_viridis_d(name="") + scale_fill_viridis_d(name="", alpha=0.25)

## ----cat_mod_boundtest, include=TRUE, fig.width=5, fig.height=4---------------
library("qgcomp")
dat2$zfactor = as.factor(dat2$z)
catfitoverall <- qgcomp.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  q = 4, family=binomial())
catfitemm <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  q = 4, family=binomial())

catfitoverallee <- qgcomp.glm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  q = 4, family=binomial())
catfitemmee <- qgcomp.emm.glm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat2,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  q = 4, family=binomial(), rr=FALSE)


anova(catfitoverall$fit, catfitemm$fit, test = "Chisq")
anova(catfitemmee, catfitoverallee)


## ----catmod_gen_letter, include=TRUE------------------------------------------
 set.seed(23)
 dat3 <- simdata_quantized_emm(
  outcometype="continuous",
  n = 300,
  corr=c(0.6, 0.5, 0.3, -0.3, -0.3, 0.0),
  b0=-2,
  mainterms=c(0.1, -0.1, 0.1, 0.0, 0.1, 0.1, 0.1),
  prodterms = c(0.2, 0.0, 0.0, 0.0, 0.2, -0.2, 0.2),
  ztype = "categorical",
  q = 4,
  yscale = 2.0                            
)

dat3$zletter = as.factor(with(dat3, ifelse(z==0, "a", ifelse(z==1, "b", "c"))))
with(dat3, table(z, zletter))

qfit_letter <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "zletter",
  q = 4, family=gaussian())

qfit_letteree <- qgcomp.emm.glm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "zletter",
  q = 4, family=gaussian())

print(qfit_letter)
print(qfit_letteree)
## output stratum specific joint effect estimate for the mixture at Zletter="a"
print(getstrateffects(qfit_letter, emmval="a"))
print(getstrateffects(qfit_letteree, emmval="a"))
## output the weights at Z=2
print(getstratweights(qfit_letter, emmval="a"))

## output predictions at Zletter="b"
print(pointwisebound(qfit_letter, emmval="b", pointwiseref = 1, alpha=0.05))
print(pointwisebound(qfit_letteree, emmval="b", pointwiseref = 1, alpha=0.05))


## ----contmod_gen, include=TRUE------------------------------------------------
 set.seed(23)
 dat3 <- simdata_quantized_emm(
  outcometype="continuous",
# sample size
  n = 100,
# correlation between x1 and x2, x3, ...
  corr=c(0.8, 0.6, 0.3, -0.3, -0.3, -0.3),
# model intercept
  b0=-2,
# linear model coefficients for x1, x2, ... at referent level of interacting variable
  mainterms=c(0.3, -0.1, 0.1, 0.0, 0.3, 0.1, 0.1),
# linear model coefficients for product terms between x1, x2, ... and interacting variable  
  prodterms = c(1.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2),
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
cor(dat3[, paste0("x", 1:7)])

## ----cont_mod_fit, include=TRUE-----------------------------------------------
qfit3 <- qgcomp.emm.glm.noboot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "CoM",
  q = 4)
qfit3
qfit3b <- qgcomp.emm.boot(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "CoM",
  B=10,
  q = 4)
qfit3b
qfit3ee <- qgcomp.emm.ee(y~x1+x2+x3+x4+x5+x6+x7,
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "CoM",
  q = 4)
qfit3ee

## ----cont_mod_bound, include=TRUE, fig.width=5, fig.height=4------------------
## output/plot the weights at CoM=0")
getstratweights(qfit3, emmval=0)
plot(qfit3, emmval=0)

## output stratum specific joint effect estimate for the mixture at CoM=0
print(getstrateffects(qfit3, emmval=0))

## output pointwise comparisons at CoM=0
print(pointwisebound(qfit3, emmval=0))


## output/plot the weights at the 80%ile of CoM
getstratweights(qfit3, emmval=quantile(dat3$CoM, .8))
plot(qfit3, emmval=quantile(dat3$CoM, .8))


## output stratum specific joint effect estimate for the mixture at the 80%ile of CoM
print(getstrateffects(qfit3, emmval=quantile(dat3$CoM, .8)))

## output pointwise comparisons at at the 80%ile of CoM
print(pointwisebound(qfit3, emmval=quantile(dat3$CoM, .8), pointwiseref = 2))


## ----cont_mod_bound_boot, include=TRUE, fig.width=5, fig.height=4-------------
## output stratum specific joint effect estimate for the mixture at CoM=0
#print(getstrateffects(qfit3b, emmval=0)) # commmented out to reduce output
print(getstrateffects(qfit3ee, emmval=0))

## output pointwise comparisons at CoM=0
#print(pointwisebound(qfit3b, emmval=0)) # commmented out to reduce output
print(pointwisebound(qfit3ee, emmval=0))

## output stratum specific joint effect estimate for the mixture at the 80%ile of CoM
#print(getstrateffects(qfit3b, emmval=quantile(dat3$CoM, .8))) # commmented out to reduce output
print(getstrateffects(qfit3ee, emmval=quantile(dat3$CoM, .8)))

## output pointwise comparisons at at the 80%ile of CoM
#print(pointwisebound(qfit3b, emmval=quantile(dat3$CoM, .8))) # commmented out to reduce output
print(pointwisebound(qfit3ee, emmval=quantile(dat3$CoM, .8)))

## plot the pointwise effects at CoM=0 (compare bootstrap with estimating equations approach)
## note that the bootstrapped version is not valid because of only using 10 bootstrap iterations
ppb <- plot(qfit3b, emmval=0, suppressprint = TRUE, geom_only = TRUE)
ppee <- plot(qfit3ee, emmval=0, suppressprint = TRUE, geom_only = TRUE)

ppb$point$data$col = "Bootstrap"
ppb$errorbar$data$col = "Bootstrap"
ppee$point$data$col = "EE"
ppee$errorbar$data$col = "EE"

ppee$point$data$x = ppee$point$data$x + 0.01
ppee$errorbar$data$x = ppee$errorbar$data$x + 0.01


ggplot() + ppb + ppee + scale_color_grey(name="Pointwise estimates, 95% CI")


## ----cont_mod_nlfit, include=TRUE---------------------------------------------
qfit3bnl <- qgcomp.emm.glm.boot(y~x1+x2+x3+x4+x5+x6+x7 + x1*(x2 + x3 + x4 + x5 + x6 +x7),
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "CoM",
  B = 10, # set higher in practice
  q = 8, degree= 2)

qfit3bnlee <- qgcomp.emm.glm.ee(y~x1+x2+x3+x4+x5+x6+x7 + x1*(x2 + x3 + x4 + x5 + x6 +x7),
  data = dat3,
  expnms = paste0("x", 1:7),
  emmvar = "CoM",
  q = 8, degree= 2)

qfit3bnl
qfit3bnlee

## ----cont_mod_nl_plot_base, include=TRUE, fig.width=6, fig.height=4-----------
print(pointwisebound(qfit3bnlee, emmval=-1))
print(pointwisebound(qfit3bnlee, emmval=-1))
print(pointwisebound(qfit3bnlee, emmval=1))
plot(qfit3bnlee, emmval=-1, modelband=TRUE, pointwiseref=4)
plot(qfit3bnlee, emmval=1, modelband=TRUE, pointwiseref=4)

## ----cont_mod_nl_plot, include=TRUE, fig.width=6, fig.height=4----------------
library(ggplot2)
plotdata = data=data.frame(q=qfit3bnlee$index, ey=qfit3bnlee$y.expected, modifier=qfit3bnlee$emmvar.msm)
ggplot() + 
         geom_point(aes(x=q, y=ey, color=modifier), data=plotdata) + 
         geom_point(aes(x=q, y=ey), color="purple", data=plotdata[plotdata$modifier>1, ], pch=1, cex=3) + 
         geom_smooth(aes(x=q, y=ey), se=FALSE, color="purple", data=plotdata[plotdata$modifier>1, ], method = 'loess', formula='y ~ x') + 
         geom_smooth(aes(x=q, y=ey), se=FALSE, color="red", data=plotdata[plotdata$modifier < -1, ], method = 'loess', formula='y ~ x') + 
         geom_point(aes(x=q, y=ey), color="red", data=plotdata[plotdata$modifier < -1, ], pch=1, cex=3) + 
  theme_classic() + 
  labs(y="Expected outcome", x="Quantile score value (0 to q-1)") + 
  scale_color_continuous(name="Value\nof\nmodifier")
         


## ----survmod_gen, include=TRUE------------------------------------------------
 set.seed(23)
 dat4 <- simdata_quantized_emm(
  outcometype="survival",
# sample size
  n = 200,
# correlation between x1 and x2, x3, ...
  corr=c(0.8, 0.6, 0.3, -0.3, -0.3, -0.3),
# model intercept
  b0=-2,
# linear model coefficients for x1, x2, ... at referent level of interacting variable
  mainterms=c(0.0, -0.1, 0.1, 0.0, 0.3, 0.1, 0.1),
# linear model coefficients for product terms between x1, x2, ... and interacting variable  
  prodterms = c(0.1, 0.0, 0.0, 0.0, -0.2, -0.2, -0.2),
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

cor(dat4[, paste0("x", 1:7)])

## ----survival, include=TRUE---------------------------------------------------
qfit4 <- qgcomp.emm.cox.noboot(survival::Surv(time, d)~x1+x2+x3+x4+x5+x6+x7,
  data = dat4,
  expnms = paste0("x", 1:7),
  emmvar = "zfactor",
  q = 4)

qfit4
plot(qfit4, emmval=0)
getstratweights(qfit4, emmval=2)
getstrateffects(qfit4, emmval=2)
pointwisebound(qfit4, emmval=1)

