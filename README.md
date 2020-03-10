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

![alt text](https://github.com/tankwin08/waveformlidar/blob/master/man/figures/r_package_graphic_abstract1.png)

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
We may want to know how good of our estimation of y compared to our real waveform lidar data
```R
res<-residuals(fit1,summary=TRUE)


## you can extract the variance of the model
VarCorr(fit1, estimate = "quantile", probs = c(0.025, 0.975))

## get the waic of the model
WAIC(fit1)


```


**How to use waveformlidar**
---
As a brief introduction, we conduct Gaussian decomposition on several waveforms to extarct useful information from the wavefomrs. In addition, we also examplify the waveforms on decompostion method with GOLD and RL methods. The detailed description of these methods has been documented in our previous reserach (https://doi.org/10.1016/j.isprsjprs.2017.04.021). 
Furthermore, we also demonstrate how to generate Hyper Point Cloud (HPC) from a sample waveform dataset at a small region.
Some examples about using these functions can be found in vignettes.

***Decomposition***
---
We used two waveforms as an example to demonstrate how to conduct Gaussian, adaptive and Weibull decompositon with differnt preprocessing settings.

```R
library(data.table)
data(return)
wf<-data.table(index=c(1:nrow(return)),return)
###################################################################################################
#decomposition two examples: one is the simple waveform and another was more complex waveform
r1<- decom(wf[1,])    ##default, with smooth was applied
r2<- decom(wf[1,], smooth = FALSE)  ##use the raw waveform

##for more complicated waveform
r3<- decom(wf[182,])
##when the waveform is mixed with too much noise, we should consider use more filtering options
r4<- decom(wf[182,],smooth=TRUE,width=5)
###we also can fit the waveform with the other models such as adaptive Gaussian and Weibull
r5<-decom.adaptive(wf[182,])
r6<-decom.weibull(wf[182,])

```
We selected one of the results to demonstrate the structure of results with differnt prepocessing settings. As you can see, the r3 gave us NULL due to there is no solution when we used the Gaussian model to decompose waveform. However, with some filtering step, we can get the solution of the decompostion as shwon in r4. But this solution is not a reasonle solution since one of the As is negative. Consequently, the index number r4[[1]] and r4[[3]] return NA and NULL to indicate that the caution should be exercised on this result. 

```
r3
NULL

r4
[[1]]
[1] NA

[[2]]
       index   Estimate Std. Error   t value     Pr(>|t|)
A1       182 228.709231  7.3797916 30.991286 1.818462e-49
A2       182 -30.882612  8.3873988 -3.682025 3.964092e-04
A3       182  81.869094  5.8023865 14.109555 2.012387e-24
u1       182  41.640000  0.4578811 90.940635 1.379831e-89
u2       182  42.130641  0.8817503 47.780695 3.086295e-65
u3       182  71.680461  0.7279330 98.471233 1.243687e-92
sigma1   182  14.612510  0.4900214 29.820145 4.161742e-48
sigma2   182   3.522174  1.1290565  3.119573 2.441069e-03
sigma3   182   8.072803  0.6615389 12.203065 1.052106e-20

[[3]]
NULL

```

***Deconvolution***
---
Compared to the decomposition, the deconvolution requires more input data and additional processing steps. Generally, we should have three kinds of data as the input for the deconvolution: the return waveform (RW), corresponding outgoing pulse (OUT) and the system impulse response (SIR). The RW and OUT are directly provided by the vendor. Ideally, the SIR is obtained through the calibration process in the lab before the waveform data are collected. In our case, the NEON provided a return impulse response (RIR) which can be assumed as a prototype SIR. This system impulse was obtained through a return pulse of single laser shot from a hard ground target with a mirror angle close to nadir. Meanwhile, NEON also provided the corresponding outgoing pulse of this return impulse response (RIR_OUT). The "true" system impulse response can be obtained by deconvolving the RIR_OUT.
In this package, we provide two options for users to deal with the system impulse response (SIR). One is directly to assume the RIR as the SIR by assigning imp = RIR. Another is to obtain the SIR through deconvolving the OUT_RIR. In the function, the "true" SIR can be achieved by assigning imp = RIR and imp_out = OUT_RIR.

```R
data(return)
data(outg)  ###corresponding outgoing pulse of return
data(imp)  ##The impulse function is generally one for the whole study area or
data(imp_out) ##corresponding outgoing pulse of imp
i=1
re<-return[i,]
out<-outg[i,]
imp<-imp
imp_out<-imp_out

### option1: to obtain the true system impluse response using the return impluse repsonse (imp) and corresponding outgoing pulse (imp_out)
gold0<-deconvolution(re = re,out = out,imp = imp,imp_out = imp_out)
rl0<-deconvolution(re = re,out = out,imp = imp,imp_out = imp_out,method = "RL")

###option2: assume the return impluse repsonse RIP is the system impulse reponse (SIR)
gold1<-deconvolution(re = re,out = out,imp = imp)
rl1<-deconvolution(re = re,out = out,imp = imp,method="RL",small_paras = c(30,2,1.5,30,2,2))
plot(gold1,type="l")
lines(rl1,col="red")
```
***Method comparison***
---
To visually compare the decompostion and deconvolution results at the individual waveform level, we made a plot of three wavefroms with these three methods (Gaussian decomposition, Gold and RL deconvolution).

****Individual waveform level****

![alt text](https://github.com/tankwin08/waveformlidar/blob/master/man/figures/README_decompostion%26deconvolution_example.png)


****Point Cloud level****

![alt text](https://github.com/tankwin08/waveformlidar/blob/master/man/figures/README_point_cloud_comparison-min.png)

You also can find these results from our previous reserach: 
Tan Zhou*, Sorin C. Popescu, Keith Krause, Ryan D. Sheridan, and Eric Putman, 2017. Gold-A novel deconvolution algorithm with optimization for waveform LiDAR processing. ISPRS Journal of Photogrammetry and Remote Sensing 129 (2017): 131-150. https://doi.org/10.1016/j.isprsjprs.2017.04.021

**Hyper point cloud**
---
The routine for extracting FW LiDAR information is to convert part of waveform signals to discrete points with the decomposition or deconvolution methods, which has been proven useful for tree species identification, forest inventory, and biomass estimation. However, most of the intensity information inherent in waveforms is being ignored with the conventional methods that undoubtedly degrades the value of FW LiDAR data. Moreover, the complicated waveform processing steps perplex users and further hinder the extensive use of FW LiDAR data for vegetation characterization. To tackle these challenges, we directly convert all raw waveform signals into points to form a point cloud, named the HPC, for subsequent analysis. A HPC is a set of data points converted from all waveform signals along the pulse path by combing geo-reference information (black) with raw waveform data (blue). 

![alt text](https://github.com/tankwin08/waveformlidar/blob/master/man/figures/hyper_point_cloud_graphic_abstract-min.jpg)

```R
data(geo)
data(return)

geo$index<-NULL
colnames(geo)[1:8]<-c("x","y","z","dx","dy","dz","or","fr")
hpc<-hyperpointcloud(waveform=return,geo=geo)
```

**Individaul tree waveform voxelization**
---
The principle behind the voxel is that the neighborhood points shared the similar characteristic and the information within the homogenous unit can be represented by one quantity or one voxel.  The following example shows how to voxelize data from the HPC. The main parameter of this function is the voxel size (res) which require you to assign a vector containing three values to represent voxel size at the X, Y and Z directions. Analogous to the waveformgrid, we also can generate the quantile intensity in each voxel by adding quan argument.
```R
voxr<-waveformvoxel(hpc,res=c(1,1,0.15))
```
Here is one simple example to conduct waveform voxelization using one individual tree.
![alt text](https://github.com/tankwin08/waveformlidar/blob/master/man/figures/individual_tree_waveformvoxel_flat_tree1_60%25_maxi_0.8_0.8_0.15_filter.png)


**Citing waveformlidar and related software**

You are welcome to use the package. If you need more help on differnt dataset or cooperation, I'd love to contribute. Developing and maintaining open source software take authors a lot of time and effort yet often underappreciated contribution to scientific progress. Thus, whenever you are
using open source software (or software in general), please make sure to cite it
appropriately so that developers get credit for their work.

When using waveformlidar, please cite one or more of the following publications:

1. Tan Zhou *, Sorin Popescu, Lonesome Malambo, Kaiguang Zhao, Keith Krause. From LiDAR waveforms to Hyper Point Clouds: a novel data product to characterize vegetation structure. Remote Sensing 2018, 10(12), 1949; https://doi.org/10.3390/rs10121949

2. Tan Zhou*, Sorin C. Popescu, A. Michelle Lawing, Marian Eriksson, Bogdan M. Strimbu, and Paul C. Bürkner. Bayesian and Classical Machine Learning Methods: A Comparison for Tree Species Classification with LiDAR Waveform Signatures. Remote Sensing 10, no. 1 (2017): 39. doi:10.3390/rs10010039

3. Tan Zhou*, and S.C. Popescu, 2017. Bayesian decomposition of full waveform LiDAR data with uncertainty analysis. Remote Sensing of Environment 200 (2017): 43-62. http://dx.doi.org/10.1016/j.rse.2017.08.012

4. Tan Zhou*, Sorin C. Popescu, Keith Krause, Ryan D. Sheridan, and Eric Putman, 2017. Gold-A novel deconvolution algorithm with optimization for waveform LiDAR processing. ISPRS Journal of Photogrammetry and Remote Sensing 129 (2017): 131-150. https://doi.org/10.1016/j.isprsjprs.2017.04.021

5. Tan Zhou*, Sorin Popescu. waveformlidar: An R Package for Waveform LiDAR Processing and Analysis. Remote Sensing. 2019, 11(21), 2552. https://doi.org/10.3390/rs11212552

**What is the best way to ask a question or propose a new feature?**
---
To propose a new feature or report a bug, please open an issue on github (https://github.com/tankwin08/waveformlidar/issues). Of course, you can always write me an email (tankchow12@gmail.com).
