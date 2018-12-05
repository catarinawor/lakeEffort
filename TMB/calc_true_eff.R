#=========================================================
# R code to run lake effort model
# version 1 - tranlated from JAGS code ImputeEff_h_1.txt
# model presented in van Poorten et al. 2015 Imputing 
# recreational angling effort from time-lapse cameras 
# using an hierarchical Bayesian model
# Author: Catarina Wor
# Date: Feb 27th 2018
#=========================================================

source("~/Documents/lake_effort/R/directories.R")

setwd(model_dir)

#source("alllakes_TMB_data_last.R")

library(tmbstan)

#library
if ("TMB" %in% row.names(installed.packages())  == FALSE) {
	install.packages("TMB") 
}
library(TMB)

setwd(TMB_dir)
#compile cpp

compile("calc_true_eff.cpp")
################################
#Make sure compilation exit with no errors -- 0
################################

dyn.load(dynlib("calc_true_eff"))

#data input -- need to process data first

load("TMB_data.Rdata")
names(TMB_data)


param <- list( md=0,pd=0.5,
			   ml=0,pl=1.5,
			    ma=0,pa=0.01,
				d=rep(0,TMB_data_jags$nlc),
				L=rep(.4, TMB_data_jags$nlc),
				a= rep(0, TMB_data_jags$nlc) 
				)

#param <- list( md=-0.2,pd=0.5,
#			   ml=1.1,pl=1.5,
#			    ma=0.6,pa=0.01,
#				d=rep(0,TMB_data$zz),
#				L=rep(.4, TMB_data$zp),
#				a= rep(0.5, TMB_data$pp) 
#				)

TMB_data$npp
TMB_data$ppc

names(TMB_data)
sort(unique(TMB_data$lc_zz))

obj <- MakeADFun(data=TMB_data_jags,
				parameters=param,
				random=c("d","L","a"),
				DLL="calc_true_eff")
 
obj$fn()
obj$gr()
opt<-nlminb(obj$par,obj$fn,obj$gr)
sdrep<-sdreport(obj)
(rep<-obj$report())

#init=c(md=0,pd=0.5,
#			   ml=0,pl=1.5,
#			    ma=0,pa=0.01,
#				d=rep(0,TMB_data_jags$nlc),
#				L=rep(.4, TMB_data_jags$nlc),
#				a= rep(0, TMB_data_jags$nlc))

fitmcmc <- tmbstan(obj, chains=1,
              iter=100000, init="random" ,
              lower=c(-5.0,0.0001,
              		 -5.0,0.0001,
              		 -5.0,0.0001,
              		 rep(-3,TMB_data_jags$nlc),
              		 rep(0, TMB_data_jags$nlc),
              		 rep(-3,TMB_data_jags$nlc)),
              upper=c(5.0,100.0,
              		 5.0,100.0,
              		 5.0,100.0,
              		 rep(20,TMB_data_jags$nlc),
              		 rep(20, TMB_data_jags$nlc),
              		 rep(20,TMB_data_jags$nlc)),
               control = list(adapt_delta = 0.98))


mc <- extract(fitmcmc, pars=names(obj$par),
              inc_warmup=FALSE, permuted=FALSE)
fit_summary <- summary(fitmcmc)

data.frame(TMB=fit_summary$summary[1:6,1],
JAGS=c(jags_output$BUGSoutput$mean$md,
jags_output$BUGSoutput$mean$pd,
jags_output$BUGSoutput$mean$ma,
jags_output$BUGSoutput$mean$pa,
jags_output$BUGSoutput$mean$ml,
jags_output$BUGSoutput$mean$pl))


r

names(rep)

names(rep$value)

unique(names(rep$value))
names(rep$par.fixed)
names(rep$par.random)

rep <- sdreport(obj)
prob <- rep$value[names(rep$value)=="prob"]
lambpp <- rep$value[names(rep$value)=="lambpp"]
lambda  <- rep$value[names(rep$value)=="lambda"]

delta <- rep$value[names(rep$value)=="delta"]
Di<-rep$value[names(rep$value)=="d"]

length(Di)



#========================================================================
#imputation


    # imputations for zero camera effort
niz<- nrow(imputation_zero)
zc_ind<-zero_cam$lakeyear_match
zimp_ind<-imputation_zero$lakeyear_match
lc_iz    = imputation_zero_lakecams

length(unique(lc_iz))
zero_cam_lakecams

match(zero_cam_lakecams,imputation_zero_lakecams,)
unique(zc_ind)
unique(zimp_ind)


lakeyr_

for(aa in 1:length(delta)){


}



for(yy in 1:niz){

		delta_iz(yy) =((1.0)/((1.0)+exp(-Di(lc_iz(yy))))) *pseen_iz(yy)
      #delta_iz(yy) = delta[something[yy]] ;
      pos(yy) = rbinom(delta_iz(yy));

      lambda_iz(yy) = exp(L(lc_iz(yy)));
      eff((yy)) = rpois(lambda_iz(i))  # use lambda to predict a positive effort given zero camera effort

      Iz[yy] = pos(yy)*eff(yy)/ncam_iz(yy)

    }



par(mfrow=c(1,3))
plot(lake_cam_zero,prob)
plot(lake_cam_falsezero,lambda)
plot(lake_cam_positive,lambpp)





save(pl,plsd,file="rw.RData")