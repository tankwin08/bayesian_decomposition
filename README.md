**Introduction**
---
The subsquent section is mainly to show how to quantify the uncertaity of waveform processing uinsg Bayesian method. The detailed description of the approach can refer to https://www.researchgate.net/publication/319018940_Bayesian_decomposition_of_full_waveform_LiDAR_data_with_uncertainty_analysis

In the domain of LiDAR applications, the observations or data are
inherently subject to various errors such as system setting, system calibration,
and range measurement errors. Additionally, the LiDAR vendors often do not clearly state what
errors are considered when the data are provided. Thus, uncertainty of
"truth" is ubiquitous and inherently present in the realities of LiDAR
data modeling.

Furthermore, the models used here are based on the non-linear functions that
generally suffer from problem of non-uniqueness, which can generate different parameter combinations given the same observational data and model, or several models can fit observational data at the expense of violating the physical meaning and theoretical assumption of the "real" model.

These problems are more evident for the sophisticated models with multiple peak components in the waveform decomposition. Thus, estimating model uncertainty is
imperative for an in-depth understanding of information derived from
data and the estimation accuracy.


Overall,the Bayesian decomposition is the intersection of waveform decompostion, Bayesian and non-linear dynamics.

![alt text](https://github.com/tankwin08/bayesian_decomposition/blob/master/man/figures/introduction1.png)

In my paper, combined these sections to conduct decompostion of waveform processing over a large region and also conduct uncertianty analysis of our processing at the parameter, derived point cloud and surface model levels.

![alt text](https://github.com/tankwin08/bayesian_decomposition/blob/master/man/figures/bayesian_uncertainty.jpg)

https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html#plotting-point-summaries-and-intervals


https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html

**How to implement it**

The following section will give you a simple example of the Bayesian decompostion using waveform lidar dataset (or any csv file) from waveformlidar package.

---
### Import packages 
if you don't have these packages, please install it using function install.packages()

```R
library(waveformlidar)
library(data.table)
library(caTools)
library(brms)
```
### Preparation 

This function will be needed for defining the fromula for non-linear later when we used brms package.

```
gennon <- function(A, u, sig) {
  n <- seq_along(A)
  fs <- paste( paste(paste0("A",n,"+","u",n,"+","sigma",n),collapse=" + "),
               "~ 1"
  )
  return(formula=as.formula(fs))
}
``` 
### Data
We used the dataset from waveformlidar package, you also can used other dataset after converting them into csv or data frame. Since waveform lidar data can be in various format.

```

data("return") ##waveform lidar data in data frame format only intensity, no other information

return1<- data.frame(id =1:nrow(return), return) ##give a index to these data, each row represent one individual waveform

```

### Model preparation

```R

####determien the number of peak for each waveform

nps<-apply(return,1,npeaks)
ite<-as.numeric(names(table(nps)))## get the unique number of peaks for each waveform or the whole dataset
il<-length(ite)
ori_priors<- apply(return1,1,peakfind)

```
### select a simple waveform to explore the Bayesian decomposition

```R
iid = 1
pn<- nps[iid]

####get the priors from all distribution for this peak
priors<-ori_priors[iid]
  
####to collect all information about these prior information
par_dat<- data.frame(do.call(rbind,priors))


##prepare the data

y<-as.numeric(return1[iid,]);
index<-y[1]
y<-y[-c(1)]
y[y==0]<-NA


y<-y-min(y,na.rm = T)+1
y<-runmean(y,3,"C")
z<-pn

df<-data.frame(x=seq_along(y),y)

## this will be used for prior and formula for the brm model
gi<-par_dat$A
gu<-par_dat$u
gsd<-par_dat$sigma

####automstically generate a formula
init<- gennls(gi, gu, gsd)
f1<-init$formula
prior<-init$start
  
```

###preprare the parameter setting up for the model

We need to assign some probability distribution to the parameters of interest. Generally the normal distribution is used to define the prior distribution of paramaters. In waveformlidar package, the peakfind can be used to roughly estimate the possible parameters for Gaussian decompostion. Here we used this function to obtain the prior information of parameters.


```R
  nn<-names(prior);l<-length(prior)

  vals<-as.numeric(prior); sp<-c()
  
  #width<-c(isd,usd,ssd) ####use the information summaried from above
  
  width<-c(rep(8,z),rep(5,z),rep(3,z))  ####specify the information based on my experience
  
  for (jj in 1:l){
    sp<-c(sp,paste("normal","(",vals[jj],",",width[jj],")",sep=""))
  }

  non1<-gennon(gi, gu, gsd)
  
  prior_final <- c(set_prior(sp[1], nlpar = nn[1]),
                 set_prior(sp[2], nlpar = nn[2]),
                 set_prior(sp[3], nlpar = nn[3]))


```
You can set different prior distribution to model the waveform such as uniform distribution of these prior information. OR for some simple waveform with no noise, you can get results estimated without assigning the prior but take longer time to obtain the results compared to the model with prior assigned. For the complicated waveform, it bettwe to use prior distributio as it may give you some unreasonable results.

## Bayesian Model building 

```R

  fit1<- brm(bf(init$formula,non1,nl=TRUE),
             data = df,
             prior = prior_final,
             chains = 1,iter = 10000,
             thin=3,
             warmup=2000,control = list(adapt_delta = 0.9))

```
## Get the parameters from our Bayesian estimation

```R
 fit1 
 sfpars<-summary(fit1)$fixed
 
 sfpars
 
                   Estimate Est.Error  l-95% CI  u-95% CI Eff.Sample      Rhat
A1_Intercept     373.02688 6.0129526 361.11689 384.43429   2535.241 1.0000382
u1_Intercept      36.59341 0.3134611  35.96919  37.19941   2659.896 1.0015244
sigma1_Intercept  10.78561 0.2978568  10.21403  11.37284   2667.000 0.9996936

```
The above will give you the estimation of parameters and their corresponding 95% credible interval (CI). This is different from frequentist statistics, we generally called it confidence interval (CI). The difference of them can refer to (https://stats.stackexchange.com/questions/2272/whats-the-difference-between-a-confidence-interval-and-a-credible-interval).

### Rhat and Eff.sample
We measured the model's convergence using the potential scale
reduction factor, named Rhat (R???), which is a statistical criterion to test
how well the Markov Chains are mixing, or moving around the parameter
space. R??? close to one indicates convergence, while high R??? value
implies that we should run a longer chain to improve convergence to
the stationary distribution. The effective sample size was also generated
to represent the equivalent number of independent iterations of the
chain. It is a criterion for the estimation efficiency. Generally, the
higher the effective sample size, the more reliable estimates can be
achieved (Gelman et al., 2014).

## Visualization of our results

The following will give you a brief overview of Bayesian decompositon results of the paramters. The left figures the distribution of estimated paramters and right ones are the MCMC sampling process. The sigma (not strat with b) is show the distribution of residual.


```R
## You can see one paramter
plot(fit1, pars = "^b",parameters="b_A1_Intercept")

##You can see all of paramters
plot(fit1)


##plot y and predicted y with uncertainty

##get the posterior samples of parameters

post_sam<-as.mcmc(fit1,pars=c("A1","u1","sigma1"))
par_samples<-data.frame(post_sam[[1]])

lr<-ncol(par_samples)/3


y1<-df[!is.na(df$y),]$y
x<-seq_along(y1)


plot(seq_along(y1),y1,type="l",col="black",xlab="(a) Time(ns)",ylab="Intensity",cex.lab=2.5,lwd=3.5,axes=F)
axis(side=1, at=seq(0,length(x)+10,15),cex.axis=2.5)
axis(side=2, at=seq(0,max(y1,na.rm=T)+5,10),cex.axis=2.5)

## add unceetainty here
for (i in 1:nrow(par_samples)){
	lines(x,par_samples[i,1] * exp(-abs(x - par_samples[i,2])^2/(2 * par_samples[i,3]^2)),lty=2,col="gray",lwd=0.2)
	}


g1<-gpars[,1];A1<-g1[1];u1<-g1[2];sigma1<-g1[3];

lines(x,A1 * exp(-abs(x - u1)^2/(2 * sigma1^2)),lty=2,col="red",lwd=3)

legend("topleft",legend=c("SW","FW"),col=c(1,2),lty=c(1,2),lwd=3,cex=2,box.lwd="white")


##add smoothed line
lines(x,y1,lwd=3,col="black")
###add text to tell the estimate error

text(max(x)-8,max(y1)-8,"WAIC= 402.8",cex=2.5)
text(max(x)-8,max(y1)-20,"SE= 8.3",cex=2.5)


```
![alt text](https://github.com/tankwin08/bayesian_decomposition/blob/master/man/figures/visualization_y_uncertainty.png)


We may want to know how good of our estimation of y compared to our real waveform lidar data
```R
res<-residuals(fit1,summary=TRUE)


## you can extract the variance of the model
VarCorr(fit1, estimate = "quantile", probs = c(0.025, 0.975))

## get the waic of the model
WAIC(fit1)


```



