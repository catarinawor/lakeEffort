################################################################################
#
#set file directory name
TMB_dir<-"~/Documents/lake_effort/TMB"

# load in libraries

library(tmbstan)
library(TMB)
library(R2jags)

#Compare Jags and TMB results
#
setwd(TMB_dir)
load("jags_data.Rdata")

# Run all possible models and find model with lowest DIC
 # Specify locations of model files
model_file <-  "/Users/catarinawor/Documents/lake_effort/JAGS/input/ImputeEff_h_1_noimp.txt"#,  # no covariates         

# Parameters for model
parsj <-  c("md","pd","ml","pl","ma","pa","d_track","L","a")
                  
# Initial parameter values for each model
initsj <-  function() {list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5 )
                          }
    # Deviance information criterion for each model
nburn     = 100000
nits      = 200000
nchain    = 2
nthin     = 2

# Run jags
jags_output <- jags( data  = jags_data,
                     inits              = initsj,
                     parameters.to.save = parsj,
                     model.file         = model_file,
                     n.chains           = nchain,
                     n.iter             = nburn+nits,
                     n.burnin           = nburn,
                     n.thin             = nthin,
                     DIC                = TRUE )
#will spit out a bunch of warnings re unused variables, just ignore them.    
#===============================================================
#TMB

setwd(TMB_dir)
#compile cpp

compile("calc_true_eff.cpp")
compile("calc_true_eff_sim.cpp")
################################
#Make sure compilation exit with no errors -- 0
################################

dyn.load(dynlib("calc_true_eff_sim"))

#data input -- need to process data first

load("TMB_data.Rdata")

#TMB_data$lc_zz<-match(TMB_data_jags$lc_zz,sort(unique(TMB_data_jags$lc_zz)))-1
#TMB_data$lc_zp<-match(TMB_data_jags$lc_zp,sort(unique(TMB_data_jags$lc_zz)))-1
#TMB_data$lc_pp<-match(TMB_data_jags$lc_pp,sort(unique(TMB_data_jags$lc_zz)))-1

#TMB_data$nlc<-length(sort(unique(TMB_data$lc_zz)))
#TMB_data$nlc<-50

names(TMB_data)

param <- list( md=0,pd=0.5,
         ml=1,pl=7.5,
          ma=0,pa=8.5,
        d=rep(0,TMB_data$nlc),
        L=rep(.4, TMB_data$nlc),
        a= rep(0, TMB_data$nlc) 
        )

obj <- MakeADFun(data=TMB_data,
        parameters=param,
        random=c("d","L","a"),
        DLL="calc_true_eff_sim")
 
obj$fn()
obj$gr()
opt<-nlminb(obj$par,obj$fn,obj$gr)
sdrep<-sdreport(obj)
(rep<-obj$report())
names(rep)
fitmcmc <- tmbstan(obj, chains=2,
              iter=100000, thin=2,
              init=list( list( md=0,pd=0.5,
                              ml=0,pl=.5,
                              ma=0,pa=0.5,
                              d=rep(0,TMB_data$nlc),
                              L=rep(.4, TMB_data$nlc),
                              a= rep(0, TMB_data$nlc) 
                        ),
                        list( md=0,pd=.01,
                              ml=0,pl=.002,
                              ma=0,pa=.001,
                              d=rep(-.1,TMB_data$nlc),
                              L=rep(.1, TMB_data$nlc),
                              a= rep(-0.1, TMB_data$nlc) 
                        )),

              #"random" ,
              lower=c(-5.0,0.0001,
                   -5.0,0.0001,
                   -5.0,0.0001,
                   rep(-3,TMB_data$nlc),
                   rep(0, TMB_data$nlc),
                   rep(-3,TMB_data$nlc)),
              upper=c(5.0,100000,
                   5.0,100000,
                   5.0,100000,
                   rep(20,TMB_data$nlc),
                   rep(20, TMB_data$nlc),
                   rep(20,TMB_data$nlc)),
               control = list(adapt_delta = 0.95))


mc <- extract(fitmcmc, pars=names(obj$par),
              inc_warmup=FALSE, permuted=FALSE)
fit_summary <- summary(fitmcmc)




#======================================
#Imputation simulation


nsims<-1000

sumEffsims<-matrix(NA,ncol=nsims,nrow=TMB_data$nallEff_ind)

for(i in 1:nsims){

  imputa<-obj$simulate(complete=TRUE)

  imputa$aggIp<-( aggregate(imputa$Ipfull,list(TMB_data$pos_agg_ind),sum))$x
  sumEff <- numeric(length=TMB_data$nallEff_ind)

  for( n in 1:TMB_data$nallEff_ind){

     if(!is.na(TMB_data$allEff_ind_pos[n])){
      sumEff[n] <- sumEff[n]+imputa$aggIp[TMB_data$allEff_ind_pos[n]+1]
     }

     if(!is.na(TMB_data$allEff_ind_zero[n])){
      sumEff[n] <- sumEff[n]+imputa$aggIz[TMB_data$allEff_ind_zero[n]+1]
     }
  }
 
  sumEffsims[,i]<-sumEff   

}


  








