

####Bayesian decomposition

library(waveformlidar)
library(data.table)
library(caTools)
library(brms)

rm(list=ls())
####to crate the non-linear list
gennon <- function(A, u, sig) {
  n <- seq_along(A)
  fs <- paste( paste(paste0("A",n,"+","u",n,"+","sigma",n),collapse=" + "),
               "~ 1"
  )
  return(formula=as.formula(fs))
}


data("return")

return1<- data.frame(id =1:nrow(return), return)




##used for collecting data
ga<-NULL
pa<-NULL

####determien the number of peak we should use

nps<-apply(return,1,npeaks)
ite<-as.numeric(names(table(nps)))## get the unique number of peaks for each waveform or the whole dataset
il<-length(ite)
ori_priors<- apply(return1,1,peakfind)

for (j in 1:il){
  pn<-ite[j]    ###to obtain 
  #############################################
  iid<-which(nps==pn) ###for just one peak, and then you do the second peaks
  
  ################################################################################
  ####get the priors from all distribution for this peak
  priors<-ori_priors[iid]
  
  ####to collect all information about these prior information
  
  par_dat<- do.call(rbind,priors)
  
  ##get an average priors for all data
  
  group<- rep(c(1:pn), (nrow(par_dat)/pn))
  par_dat<- data.table(data.frame(par_dat,group = group))
  colnames(par_dat)[1:4]<- c("index","A","u","sig")
  
  
  ##get the boundary of parameters
  dd = par_dat[,list(a_quan = quantile(A, c(0.05,0.5,0.95)),
                u_quan = quantile(u,c(0.05,0.5,0.95)),
                sig_quan = quantile(sig,c(0.05,0.5,0.95)),
                a_sd = sd(A),
                u_sd = sd(u),
                sig_sd =sd(sig)
                ),by=group]
  

  
  ####for parameters limit, we can define here
  ####based on our experiment, this step is very important, it can influence our result a lot,
  ####we will summarized the data and determine the prior's ub and lb
  lbound_ind<- seq(1,nrow(dd),3)
  ubound_ind<- seq(3,nrow(dd),3)
  
  gil<-dd[lbound_ind,]$a_quan;giu<-dd[ubound_ind,]$a_quan
  
  gul<-dd[lbound_ind,]$u_quan;guu<-dd[ubound_ind,]$u_quan
  
  gsdl<-dd[lbound_ind,]$sig_quan;gsdu<-dd[ubound_ind,]$sig_quan
  
  ##################################################################################
  
  
  nj<-length(iid)
  ##############to create the first object
  id0<-as.numeric(iid[1]);
  y<-as.numeric(return1[id0,]);
  index<-y[1]
  #y0<-dat[id0,];y0[y0==0]<-NA
  #y0<-as.numeric(unlist(y0));index<-y0[1]
  y<-y[-c(1)]
  y[y==0]<-NA
  ###when for direct decomposition
  y<-y-min(y,na.rm = T)+1
  y<-runmean(y,3,"C")
  z<-pn
  
  ##plot
  #plot(y,type="l")
  #points(realind,newpeak)
  
  ##################################initilize parameters
  sel_row<- seq(2,3*pn,3)
  
  gi<-dd[sel_row,]$a_quan
  gu<-dd[sel_row,]$u_quan
  gsd<-dd[sel_row,]$sig_quan
  
  ####automstically generate a formula
  init<- gennls(gi, gu, gsd)
  f1<-init$formula
  prior<-init$start
  
  df<-data.frame(x=seq_along(y),y)
  
  
  ###preprare the parameter setting up
  nn<-names(prior);l<-length(prior)

  vals<-as.numeric(prior); sp<-c()
  #width<-c(isd,usd,ssd) ####use the information summaried from above
  width<-c(rep(8,z),rep(5,z),rep(3,z))  ####specify the information based on my experience
  for (jj in 1:l){
    sp<-c(sp,paste("normal","(",vals[jj],",",width[jj],")",sep=""))
  }
  
  
  ###to create the prior parameters series

  non1<-gennon(gi, gu, gsd)
  ln<-l/3
  ####to determine which prior we should use
  if (ln==1){
    prior_final <- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                 set_prior(sp[2], nlpar = nn[2],lb=gul[1],ub=guu[1]),
                 set_prior(sp[3], nlpar = nn[3],lb=gsdl[1],ub=gsdu[1]))
    
  } else if (ln==2) {
    prior_final<-c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                   set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                   set_prior(sp[3], nlpar = nn[3],lb=gul[1],ub=guu[1]),
                   set_prior(sp[4], nlpar = nn[4],lb=gul[2],ub=guu[2]),
                   set_prior(sp[5], nlpar = nn[5],lb=gsdl[1],ub=gsdu[1]),
                   set_prior(sp[6], nlpar = nn[6],lb=gsdl[2],ub=gsdu[2]))
    
  } else if (ln==3) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                               set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                               set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                               set_prior(sp[4], nlpar = nn[4],lb=gul[1],ub=guu[1]),
                               set_prior(sp[5], nlpar = nn[5],lb=gul[2],ub=guu[2]),
                               set_prior(sp[6], nlpar = nn[6],lb=gul[3],ub=guu[3]),
                               set_prior(sp[7], nlpar = nn[7],lb=gsdl[1],ub=gsdu[1]),
                               set_prior(sp[8], nlpar = nn[8],lb=gsdl[2],ub=gsdu[2]),
                               set_prior(sp[9], nlpar = nn[9],lb=gsdl[3],ub=gsdu[3]))
  } else if (ln==4) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                               set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                               set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                               set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                               set_prior(sp[5], nlpar = nn[5],lb=gul[1],ub=guu[1]),
                               set_prior(sp[6], nlpar = nn[6],lb=gul[2],ub=guu[2]),
                               set_prior(sp[7], nlpar = nn[7],lb=gul[3],ub=guu[3]),
                               set_prior(sp[8], nlpar = nn[8],lb=gul[4],ub=guu[4]),
                               set_prior(sp[9], nlpar = nn[9],lb=gsdl[1],ub=gsdu[1]),
                               set_prior(sp[10], nlpar = nn[10],lb=gsdl[2],ub=gsdu[2]),
                               set_prior(sp[11], nlpar = nn[11],lb=gsdl[3],ub=gsdu[3]),
                               set_prior(sp[12], nlpar = nn[12],lb=gsdl[4],ub=gsdu[4]))
  } else if (ln==5) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                              set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                              set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                              set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                              set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                              set_prior(sp[6], nlpar = nn[6],lb=gul[1],ub=guu[1]),
                              set_prior(sp[7], nlpar = nn[7],lb=gul[2],ub=guu[2]),
                              set_prior(sp[8], nlpar = nn[8],lb=gul[3],ub=guu[3]),
                              set_prior(sp[9], nlpar = nn[9],lb=gul[4],ub=guu[4]),
                              set_prior(sp[10], nlpar = nn[10],lb=gul[5],ub=guu[5]),
                              set_prior(sp[11], nlpar = nn[11],lb=gsdl[1],ub=gsdu[1]),
                              set_prior(sp[12], nlpar = nn[12],lb=gsdl[2],ub=gsdu[2]),
                              set_prior(sp[13], nlpar = nn[13],lb=gsdl[3],ub=gsdu[3]),
                              set_prior(sp[14], nlpar = nn[14],lb=gsdl[4],ub=gsdu[4]),
                              set_prior(sp[15], nlpar = nn[15],lb=gsdl[5],ub=gsdu[5]))
  } else if (ln==6) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                                set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                                set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                                set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                                set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                                set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                                set_prior(sp[7], nlpar = nn[7],lb=gul[1],ub=guu[1]),
                                set_prior(sp[8], nlpar = nn[8],lb=gul[2],ub=guu[2]),
                                set_prior(sp[9], nlpar = nn[9],lb=gul[3],ub=guu[3]),
                                set_prior(sp[10], nlpar = nn[10],lb=gul[4],ub=guu[4]),
                                set_prior(sp[11], nlpar = nn[11],lb=gul[5],ub=guu[5]),
                                set_prior(sp[12], nlpar = nn[12],lb=gul[6],ub=guu[6]),
                                set_prior(sp[13], nlpar = nn[13],lb=gsdl[1],ub=gsdu[1]),
                                set_prior(sp[14], nlpar = nn[14],lb=gsdl[2],ub=gsdu[2]),
                                set_prior(sp[15], nlpar = nn[15],lb=gsdl[3],ub=gsdu[3]),
                                set_prior(sp[16], nlpar = nn[16],lb=gsdl[4],ub=gsdu[4]),
                                set_prior(sp[17], nlpar = nn[17],lb=gsdl[5],ub=gsdu[5]),
                                set_prior(sp[18], nlpar = nn[18],lb=gsdl[6],ub=gsdu[6]))
  } else if (ln==7) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                                set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                                set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                                set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                                set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                                set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                                set_prior(sp[7], nlpar = nn[7],lb=gil[7],ub=giu[7]),
                                set_prior(sp[8], nlpar = nn[8],lb=gul[1],ub=guu[1]),
                                set_prior(sp[9], nlpar = nn[9],lb=gul[2],ub=guu[2]),
                                set_prior(sp[10], nlpar = nn[10],lb=gul[3],ub=guu[3]),
                                set_prior(sp[11], nlpar = nn[11],lb=gul[4],ub=guu[4]),
                                set_prior(sp[12], nlpar = nn[12],lb=gul[5],ub=guu[5]),
                                set_prior(sp[13], nlpar = nn[13],lb=gul[6],ub=guu[6]),
                                set_prior(sp[14], nlpar = nn[14],lb=gul[7],ub=guu[7]),
                                set_prior(sp[15], nlpar = nn[15],lb=gsdl[1],ub=gsdu[1]),
                                set_prior(sp[16], nlpar = nn[16],lb=gsdl[2],ub=gsdu[2]),
                                set_prior(sp[17], nlpar = nn[17],lb=gsdl[3],ub=gsdu[3]),
                                set_prior(sp[18], nlpar = nn[18],lb=gsdl[4],ub=gsdu[4]),
                                set_prior(sp[19], nlpar = nn[19],lb=gsdl[5],ub=gsdu[5]),
                                set_prior(sp[20], nlpar = nn[20],lb=gsdl[6],ub=gsdu[6]),
                                set_prior(sp[21], nlpar = nn[21],lb=gsdl[7],ub=gsdu[7]))
  } else if (ln==8) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                                set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                                set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                                set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                                set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                                set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                                set_prior(sp[7], nlpar = nn[7],lb=gil[7],ub=giu[7]),
                                set_prior(sp[8], nlpar = nn[8],lb=gil[8],ub=giu[8]),
                                set_prior(sp[9], nlpar = nn[9],lb=gul[1],ub=guu[1]),
                                set_prior(sp[10], nlpar = nn[10],lb=gul[2],ub=guu[2]),
                                set_prior(sp[11], nlpar = nn[11],lb=gul[3],ub=guu[3]),
                                set_prior(sp[12], nlpar = nn[12],lb=gul[4],ub=guu[4]),
                                set_prior(sp[13], nlpar = nn[13],lb=gul[5],ub=guu[5]),
                                set_prior(sp[14], nlpar = nn[14],lb=gul[6],ub=guu[6]),
                                set_prior(sp[15], nlpar = nn[15],lb=gul[7],ub=guu[7]),
                                set_prior(sp[16], nlpar = nn[16],lb=gul[8],ub=guu[8]),
                                set_prior(sp[17], nlpar = nn[17],lb=gsdl[1],ub=gsdu[1]),
                                set_prior(sp[18], nlpar = nn[18],lb=gsdl[2],ub=gsdu[2]),
                                set_prior(sp[19], nlpar = nn[19],lb=gsdl[3],ub=gsdu[3]),
                                set_prior(sp[20], nlpar = nn[20],lb=gsdl[4],ub=gsdu[4]),
                                set_prior(sp[21], nlpar = nn[21],lb=gsdl[5],ub=gsdu[5]),
                                set_prior(sp[22], nlpar = nn[22],lb=gsdl[6],ub=gsdu[6]),
                                set_prior(sp[23], nlpar = nn[23],lb=gsdl[7],ub=gsdu[7]),
                                set_prior(sp[24], nlpar = nn[24],lb=gsdl[8],ub=gsdu[8]))
    
  } else if (ln==9) {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                                set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                                set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                                set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                                set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                                set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                                set_prior(sp[7], nlpar = nn[7],lb=gil[7],ub=giu[7]),
                                set_prior(sp[8], nlpar = nn[8],lb=gil[8],ub=giu[8]),
                                set_prior(sp[9], nlpar = nn[9],lb=gil[9],ub=giu[9]),
                                set_prior(sp[10], nlpar = nn[10],lb=gul[1],ub=guu[1]),
                                set_prior(sp[11], nlpar = nn[11],lb=gul[2],ub=guu[2]),
                                set_prior(sp[12], nlpar = nn[12],lb=gul[3],ub=guu[3]),
                                set_prior(sp[13], nlpar = nn[13],lb=gul[4],ub=guu[4]),
                                set_prior(sp[14], nlpar = nn[14],lb=gul[5],ub=guu[5]),
                                set_prior(sp[15], nlpar = nn[15],lb=gul[6],ub=guu[6]),
                                set_prior(sp[16], nlpar = nn[16],lb=gul[7],ub=guu[7]),
                                set_prior(sp[17], nlpar = nn[17],lb=gul[8],ub=guu[8]),
                                set_prior(sp[18], nlpar = nn[18],lb=gul[9],ub=guu[9]),
                                set_prior(sp[19], nlpar = nn[19],lb=gsdl[1],ub=gsdu[1]),
                                set_prior(sp[20], nlpar = nn[20],lb=gsdl[2],ub=gsdu[2]),
                                set_prior(sp[21], nlpar = nn[21],lb=gsdl[3],ub=gsdu[3]),
                                set_prior(sp[22], nlpar = nn[22],lb=gsdl[4],ub=gsdu[4]),
                                set_prior(sp[23], nlpar = nn[23],lb=gsdl[5],ub=gsdu[5]),
                                set_prior(sp[24], nlpar = nn[24],lb=gsdl[6],ub=gsdu[6]),
                                set_prior(sp[25], nlpar = nn[25],lb=gsdl[7],ub=gsdu[7]),
                                set_prior(sp[26], nlpar = nn[26],lb=gsdl[8],ub=gsdu[8]),
                                set_prior(sp[27], nlpar = nn[27],lb=gsdl[9],ub=gsdu[9]))
    
  } else if (ln==10) {
    prior_final<-c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                   set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                   set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                   set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                   set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                   set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                   set_prior(sp[7], nlpar = nn[7],lb=gil[7],ub=giu[7]),
                   set_prior(sp[8], nlpar = nn[8],lb=gil[8],ub=giu[8]),
                   set_prior(sp[9], nlpar = nn[9],lb=gil[9],ub=giu[9]),
                   set_prior(sp[10], nlpar = nn[10],lb=gil[10],ub=giu[10]),
                   set_prior(sp[11], nlpar = nn[11],lb=gul[1],ub=guu[1]),
                   set_prior(sp[12], nlpar = nn[12],lb=gul[2],ub=guu[2]),
                   set_prior(sp[13], nlpar = nn[13],lb=gul[3],ub=guu[3]),
                   set_prior(sp[14], nlpar = nn[14],lb=gul[4],ub=guu[4]),
                   set_prior(sp[15], nlpar = nn[15],lb=gul[5],ub=guu[5]),
                   set_prior(sp[16], nlpar = nn[16],lb=gul[6],ub=guu[6]),
                   set_prior(sp[17], nlpar = nn[17],lb=gul[7],ub=guu[7]),
                   set_prior(sp[18], nlpar = nn[18],lb=gul[8],ub=guu[8]),
                   set_prior(sp[19], nlpar = nn[19],lb=gul[9],ub=guu[9]),
                   set_prior(sp[20], nlpar = nn[20],lb=gul[10],ub=guu[10]),
                   set_prior(sp[21], nlpar = nn[21],lb=gsdl[1],ub=gsdu[1]),
                   set_prior(sp[22], nlpar = nn[22],lb=gsdl[2],ub=gsdu[2]),
                   set_prior(sp[23], nlpar = nn[23],lb=gsdl[3],ub=gsdu[3]),
                   set_prior(sp[24], nlpar = nn[24],lb=gsdl[4],ub=gsdu[4]),
                   set_prior(sp[25], nlpar = nn[25],lb=gsdl[5],ub=gsdu[5]),
                   set_prior(sp[26], nlpar = nn[26],lb=gsdl[6],ub=gsdu[6]),
                   set_prior(sp[27], nlpar = nn[27],lb=gsdl[7],ub=gsdu[7]),
                   set_prior(sp[28], nlpar = nn[28],lb=gsdl[8],ub=gsdu[8]),
                   set_prior(sp[29], nlpar = nn[29],lb=gsdl[9],ub=gsdu[9]),
                   set_prior(sp[30], nlpar = nn[30],lb=gsdl[10],ub=gsdu[10]))
    
  } else {
    prior_final<- c(set_prior(sp[1], nlpar = nn[1],lb=gil[1],ub=giu[1]),
                                 set_prior(sp[2], nlpar = nn[2],lb=gil[2],ub=giu[2]),
                                 set_prior(sp[3], nlpar = nn[3],lb=gil[3],ub=giu[3]),
                                 set_prior(sp[4], nlpar = nn[4],lb=gil[4],ub=giu[4]),
                                 set_prior(sp[5], nlpar = nn[5],lb=gil[5],ub=giu[5]),
                                 set_prior(sp[6], nlpar = nn[6],lb=gil[6],ub=giu[6]),
                                 set_prior(sp[7], nlpar = nn[7],lb=gil[7],ub=giu[7]),
                                 set_prior(sp[8], nlpar = nn[8],lb=gil[8],ub=giu[8]),
                                 set_prior(sp[9], nlpar = nn[9],lb=gil[9],ub=giu[9]),
                                 set_prior(sp[10], nlpar = nn[10],lb=gil[10],ub=giu[10]),
                                 set_prior(sp[11], nlpar = nn[11],lb=gil[11],ub=giu[11]),
                                 set_prior(sp[12], nlpar = nn[12],lb=gul[1],ub=guu[1]),
                                 set_prior(sp[13], nlpar = nn[13],lb=gul[2],ub=guu[2]),
                                 set_prior(sp[14], nlpar = nn[14],lb=gul[3],ub=guu[3]),
                                 set_prior(sp[15], nlpar = nn[15],lb=gul[4],ub=guu[4]),
                                 set_prior(sp[16], nlpar = nn[16],lb=gul[5],ub=guu[5]),
                                 set_prior(sp[17], nlpar = nn[17],lb=gul[6],ub=guu[6]),
                                 set_prior(sp[18], nlpar = nn[18],lb=gul[7],ub=guu[7]),
                                 set_prior(sp[19], nlpar = nn[19],lb=gul[8],ub=guu[8]),
                                 set_prior(sp[20], nlpar = nn[20],lb=gul[9],ub=guu[9]),
                                 set_prior(sp[21], nlpar = nn[21],lb=gul[10],ub=guu[10]),
                                 set_prior(sp[22], nlpar = nn[22],lb=gul[11],ub=guu[11]),
                                 set_prior(sp[23], nlpar = nn[23],lb=gsdl[1],ub=gsdu[1]),
                                 set_prior(sp[24], nlpar = nn[24],lb=gsdl[2],ub=gsdu[2]),
                                 set_prior(sp[25], nlpar = nn[25],lb=gsdl[3],ub=gsdu[3]),
                                 set_prior(sp[26], nlpar = nn[26],lb=gsdl[4],ub=gsdu[4]),
                                 set_prior(sp[27], nlpar = nn[27],lb=gsdl[5],ub=gsdu[5]),
                                 set_prior(sp[28], nlpar = nn[28],lb=gsdl[6],ub=gsdu[6]),
                                 set_prior(sp[29], nlpar = nn[29],lb=gsdl[7],ub=gsdu[7]),
                                 set_prior(sp[30], nlpar = nn[30],lb=gsdl[8],ub=gsdu[8]),
                                 set_prior(sp[31], nlpar = nn[31],lb=gsdl[9],ub=gsdu[9]),
                                 set_prior(sp[32], nlpar = nn[32],lb=gsdl[10],ub=gsdu[10]),
                                 set_prior(sp[33], nlpar = nn[33],lb=gsdl[11],ub=gsdu[11]))
  }
  #####how to determine the prior automatically,
  ###to implent differnt itrations and warmup
  
  if (ln==1){
    itrn<-6000;wn<-2000
  } else if (ln==2){
    itrn<-10000;wn<-4000
  } else if (ln==3){
    itrn<-14000;wn<-6000
  } else if (ln==4){
    itrn<-18000;wn<-6500
  } else if (ln==5){
    itrn<-20000;wn<-7000
  } else if (ln==6){
    itrn<-23000;wn<-7500
  } else if (ln==7){
    itrn<-25000;wn<-8000
  } else {
    itrn<-28000
  }
  
  fit1<- brm(bf(init$formula,non1,nl=TRUE),
             data = df,
             prior = prior_final,
             chains = 1,iter = itrn,
             thin=3,
             warmup=wn,control = list(adapt_delta = 0.9))
  
 
  sfpars<-summary(fit1)$fixed
  ga<-rbind(ga,cbind(index,sfpars))
  for (iij in 1:z){
    si<-1+(iij-1)*3;ei<-iij*3;sdi<-2+(iij-1)*3
    pm<-sfpars[si:ei,1];upars<-sfpars[sdi,2:6]
    pmc<-c(index,pm,upars)
    pa<-rbind(pa,pmc)
  }
  
  
  
  nj<-length(iid)
  ####20 has 4 peaks
  #sj<-round(nj*0.8)+1
  if (nj>1){
    for (i in 2:nj){
    id0<-as.numeric(iid[i]);
    y<-as.numeric(return1[id0,]);
    index<-y[1]
    y<-y[-c(1)]
    y[y==0]<-NA
    ###when for direct decomposition
    y<-y-min(y,na.rm = T)+1;
    y<-runmean(y,3,"C")
    df<-data.frame(x=seq_along(y),y)
    
    fit2<-update(fit1,newdata=df,
                 iter=itrn,warmup=wn,control = list(adapt_delta = 0.9))  ####this is very important
    sfpars<-summary(fit2)$fixed
    ###to get all paramters
    ga<-rbind(ga,cbind(index,sfpars))
    
    ###directly get some parameters
    ##make a matrix, here we will get the estimate of index, A1,u1,sigma1,u1 estimate error, u1 low95% and u1 up95%
    for (iij in 1:z){
      si<-1+(iij-1)*3;ei<-iij*3;sdi<-2+(iij-1)*3
      pm<-sfpars[si:ei,1];upars<-sfpars[sdi,c(2:6)]
      pmc<-c(index,pm,upars)
      pa<-rbind(pa,pmc)
    }
    
  }
  }
}


