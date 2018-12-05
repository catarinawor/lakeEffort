################################################################################
# imputeEffort.R
# Steven Rossi (spr1@sfu.ca)
# July 3, 2015
#Modified by :Catarina Wor for exploratory analysis
#
# Usage: imputeEffort()
#
# Input files: (a) "Camera Effort2.data"
#              (b) "Lake Info.data"
#
# Output files: (a) "MCMC.i.txt" (for i in 1:8, if run_all=TRUE)
#               (b) "Imputation_MCMC.dat"
#
# Required JAGS files: (a) "ImputeEff_h_i.txt"  (for i in 1:8, if run_all=TRUE)
#                      (b) "ImputeEff_h_impute.txt"
#
# Note: This script makes use of the dplyr R package. dplyr is built around 6 
# main functions that manipulate data frames:
#   1) mutate(), which defines new columns (often as a function of other cols)
#   2) filter(), which returns rows matching one or more conditions
#   3) select(), which returns only the columns passed as arguments
#   4) arrange(), which sorts rows in descending order for a given variable
#   5) group_by(), which groups the data frame by given variables. Subsequent
#                  operations will be performed "by group"
#   6) summarise(), which combines multiple values into a single value
# dplyr also makes use of the pipe operator (%>%) from the magrittr R package.
# Pipes allow commands to be chained by forwarding the value or result of an
# expression into the next function call/expression. For example:
#
#   b <- a %>% f1() %>% f2() %>% f3();
#
#    is equivalent to
#
#   b <- f3( f2( f1(a) ) );
#
# Enter install.packages("dplyr") at the console to install both packages.
# For more information:
#   dplyr:
#    https://cran.r-project.org/web/packages/dplyr/index.html
#    https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html
#   magittr and pipes:
#    https://cran.r-project.org/web/packages/magrittr/index.html
#    http://seananderson.ca/2014/09/13/dplyr-intro.html
#
#
################################################################################
# Load necessary packages
library(dplyr)
library(R2jags)

#this line must be changed if ran in another machine
source("~/Documents/lake_effort/R/directories.R")

setwd(input_data_dir)



#==================
#data and descripttion
#see data description requirement in the readme file
myeffdata<-read.table( "Camera Effort2.data", header=TRUE ) %>%
                as.data.frame()
head(myeffdata)

lake_dat <- read.table( "Lake Info.data", header=TRUE ) %>%
              as.data.frame() 
summary(lake_dat)





                  nchain    = 2      #,     # Number of MCMC chains
                  nburn     = 10000  #, # Number of MCMC iter's to discard
                  nits      = 20000  #, # Number of MCMC iterations
                  nthin     = 2      #,     # MCMC thinning factor
                  detail    = FALSE  #, # TRUE=daily effort imputed for one
                                     # lake. FALSE=annual effort imputed
                                     # across lakes
                  one_ly    = 6      #,     # If detail=TRUE, which lake?
                  eval_fit  = FALSE  #, # TRUE=fit to known data and
                                     # evaluate bias and precision
                  psamp     = 0.8    #,   # If eval_fit=TRUE, what proportion
                                     # of dataset to use for training?
                  n_compl   = 2      #,     # Number of lakes to use when
                                     # filling data holes
                  cutoff    = 0.3    #,   # Minimum proportion of days where
                                     # camera was working
                  mac       = 1      #,     # Mac (=TRUE) or PC (=FALSE)
                  run_all   = FALSE  #, # Run all models and check DIC?
                  imp_model = 5      # If run_all=FALSE, which model to
                                     # use for imputation?
                  

 
  #----------------------------------------------------------------------------
  #  Build effort and lake datasets
  #----------------------------------------------------------------------------

  # Build camera effort dataset
  effort_dat <- read.table( "Camera Effort2.data", header=TRUE ) %>%
                as.data.frame() %>%
                # Remove duplicated rows
                filter( !duplicated(.) ) %>%
                # Rename the columns to keep a consistent style
                rename( year         = Year,
                        name         = Name,      
                        region       = Region,
                        method       = Method,
                        hour         = Hour,
                        angler_count = TotalAnglerCount,
                        julian_day   = DoY,
                        time_of_day  = ToD ) %>%
                # Delete unneeded columns
                select( -Date, -Day, -Month, -WBID, -Time, -Include, -DoW )

  # Build lake dataset
  lake_dat <- read.table( "Lake Info.data", header=TRUE ) %>%
              as.data.frame() %>%
              # Rename the columns to keep a consistent style
              rename( area     = SA,
                      access   = Access,
                      dist_kam = DistKam,
                      dist_van = DistVan,
                      p_seen   = Pseen,
                      launch   = Launch,
                      camp     = Camp,
                      resort   = Resort ) %>%
              # Delete unneeded columns
              select( -Lake, -WBID, -Region, -Elevation, -Perim, -Mdepth,
                      -MuDepth, -Distgravel )

  # Get the individual lake names from the effort data set
  lake_names <- effort_dat %>%
                group_by( name ) %>%  # Group dataset by lake name
                summarise() %>%       # Return unique lake names alphabetically
                .$name %>%            # Return the lake name column as a vector
                as.character()        # Convert from factor to character

  # Identify lakes in effort dataset by numeric ID codes
  effort_dat$lake_ID <- match( effort_dat$name, lake_names )

  # Identify lake information for each observation in the effort dataset
  lake_dat_obs <- lake_dat[ effort_dat$lake_ID, ]

  # Merge the effort dataset with the corresponding lake information, then
  # delete unwanted/invalid observations
  effort_dat <- cbind(effort_dat,lake_dat_obs) %>%
                as.data.frame() %>% 
                # Delete Chief Gray and Roche observations. Chief Gray was 
                # posted on a trail with a motion-detector, so invalid for this
                # analysis. Roche is too large to justify without using
                # independent counts.
                filter( name!="CHIEFGRAY", 
                        name!="ROCHE" ) %>% 
                # Assign each observation a "lakeyear" code that uniquely 
                # idenfies the lake and year of the observation, as well as a
                # "datetimelake" code that uniquely identifies the date, time, 
                # and lake of the observation.
                mutate( lakeyear     = 100*year + lake_ID,
                        datetimelake = 100000*round(datetime,2) + lake_ID ) %>%
                # Delete columns that are no longer needed
                select( -name,
                        -datetime ) 
  
  # Remove objects that are no longer needed
  rm( lake_names,
      lake_dat,
      lake_dat_obs )

  #----------------------------------------------------------------------------
  #  Partition effort dataset by method of observation
  #----------------------------------------------------------------------------

  # Ground observations
  effort_ground <- effort_dat %>%
                   filter( method=="GR" ) %>%
                   rename( ground_count = angler_count ) %>%
                   select( -method )

  # Camera observations
  effort_cam <- effort_dat %>%
                filter( method=="CAM" ) %>%
                mutate( lake_camera = 10*lake_ID+location) %>%
                rename( cam_count = angler_count ) %>%
                select( -method )


  #----------------------------------------------------------------------------
  #  Create new effort dataset with ground and cam observations on same line
  #----------------------------------------------------------------------------

  # Create new dataset to hold both camera and ground observations
  # Start off with the same column names as the original effort dataset
  effort <- effort_dat[0, ]

  # Get unique locations from effort dataset
  locations <- effort_dat %>%
               group_by( location ) %>%   # Group dataset by location
               summarise() %>%            # Get unique locations
               .$location                 # Get locations as vector

  # Loop over all locations. At each location, combine the camera observations
  # with their corresponding ground observations.
  for( i in locations ) {

    # Camera observations at location i
    effort_cam_i <- filter( effort_cam, location==i )
    
    # Find the ground observations that correspond to camera observations i
    cam_index <- match( effort_ground$datetimelake, effort_cam_i$datetimelake )

    # Combine ground and camera observations
    effort_i <- cbind( effort_ground,
                       cam_count   = effort_cam_i$cam_count[cam_index],
                       lake_camera = effort_cam_i$lake_camera[cam_index] )

    # Add observations at location i to effort dataset
    effort <- rbind( effort, effort_i )
    
  }

  # Convert effort matrix to dataframe and delete datetimelake column
  effort <- effort %>%
            as.data.frame() %>%
            select( -datetimelake )


  

  #if( eval_fit ) {
  #  unif_samples <- runif( n=nrow(effort), min=0, max=1 )
  #  unif_subset  <- unif_samples < psamp
  #  valid        <- effort[!unif_subset, ] %>%
  #                  na.omit() %>%
  #                  as.data.frame()
  #  effort       <- effort[unif_subset, ]
  #}
  
  #----------------------------------------------------------------------------
  #  Partition effort dataset by type of observations (Positive, zero, etc...)
  #----------------------------------------------------------------------------

  # Positive camera observations
  pos_cam <- filter( effort, cam_count>0 )
  
  # Zero camera observations
  zero_cam <- effort %>%
              group_by( p_seen, area, weekend, time_of_day, hour, year,
                        lake_camera, lake_ID ) %>%
              arrange( p_seen, area, weekend, time_of_day, hour, year, 
                        lake_camera, lake_ID ) %>%
              filter( cam_count==0 ) %>%
              mutate( ground_count=ifelse( test = ground_count>0,
                                           yes  = 1,
                                           no   = 0 ) ) %>%
              summarise( n_zero_cam = sum(cam_count==0),
                         n_zero_eff = sum(ground_count==0) ) 
  
  # False zero: zero camera effort and positive total effort
  false_zero <- effort %>%
                filter( ground_count > 0,
                        cam_count   == 0 )
  
  # Imputation dataset: Positive camera effort but no ground counts
  imputation_dat <- effort_cam %>%
                    filter( year >= 2009 ) %>%
                    mutate( ground_count = NA ) %>%
                    arrange( lake_ID )

  summary(imputation_dat)
  
  # Total number of days in imputation dataset
  n_days_total <- imputation_dat %>%
                  group_by( julian_day ) %>%
                  summarise() %>%
                  nrow()

  # Unique lakeyears
  lakeyears <- effort_cam %>%
               group_by( lakeyear ) %>%
               summarise() %>%
               .$lakeyear

  # Number of days per lakeyear
  for( i in lakeyears ) {
    n_days_i <- imputation_dat %>%
                filter( lakeyear == i ) %>%
                group_by( julian_day ) %>%
                summarise() %>%
                nrow()
    
    # Proportion of days where camera was working
    days_prop <- n_days_i / n_days_total

    # If proportion of days where camera was working is too small, exclude
    # those observations from the imputation dataset
    if( days_prop < cutoff ) {
      imputation_dat <- filter( imputation_dat, lakeyear!=i )

      lakeyears <- lakeyears[ lakeyears!=i ]
    }

  }
  
  # Remove original effort dataframe
  rm(effort)

  if( eval_fit )
    imputation_dat <- valid

  lakeyears <- sort(unique(imputation_dat$lakeyear))

  # Number of lakeyears
  n_lakeyears <- length( lakeyears )
  

  #----------------------------------------------------------------------------
  #  Partition imputation dataset and record metadata
  #----------------------------------------------------------------------------

  # Sort imputation dataset by lakeyear
  imputation_dat <- arrange( imputation_dat, lakeyear )
  
  # Update lakeyear tags to range from 1 to n_lakeyears
  imputation_dat$lakeyear <- match(imputation_dat$lakeyear,lakeyears)
  false_zero$lakeyear     <- match(false_zero$lakeyear,lakeyears)
  pos_cam$lakeyear        <- match(pos_cam$lakeyear,lakeyears)

  # Create codes to differentiate observations
  imputation_dat <- imputation_dat %>%
                    mutate( lakeyearday     = 1000*lakeyear +
                                                   julian_day,
                            dayhourlakeyear = 100000*julian_day +
                                                1000*hour +
                                                     lakeyear )

  # If detail is TRUE then daily effort is imputed for a single lake
  # If detail is FALSE then annual effort is imputed across all lakes
  if( detail ) {
    imputation_dat <- filter(imputation_dat,lakeyear==one_ly)
    false_zero     <- filter(false_zero,    lakeyear==one_ly)
    pos_cam        <- filter(pos_cam,       lakeyear==one_ly)
  }

  # Zero camera counts in imputation dataset
  imputation_zero <- imputation_dat %>%
                     filter( cam_count==0 ) %>%
                     group_by( dayhourlakeyear ) %>%
                     mutate( n_cam_per_dhly = length(dayhourlakeyear)) %>%
                     group_by( lakeyearday ) %>%
                     mutate( day_index = 1:length(lakeyearday) )
  
  # Positive camera counts in imputation dataset
  imputation_pos <- imputation_dat %>%
                    filter( cam_count > 0 ) %>%
                    group_by( dayhourlakeyear ) %>%
                    mutate( n_cam_per_dhly = length(dayhourlakeyear)) %>%
                    group_by( lakeyearday ) %>%
                    mutate( day_index = 1:length(lakeyearday) )

  # Number of days to impute for
  n_days <- length(unique(imputation_dat$julian_day))
  
  # Matrix to store the days for imputation
  all_days <- matrix( data=1, nrow=n_lakeyears, ncol=n_days )

  # Vector to store the number of days to impute for each lake-year
  n_days_per_lakeyear <- numeric( length=n_lakeyears )

  # Final day to impute for
  max_julian_day <- max( imputation_dat$julian_day )

  # Matrices to store metadata
  pos_days       <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears )
  pos_nobs       <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears )
  pos_lakeyears  <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears )
  zero_nobs      <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears ) 
  zero_lakeyears <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears )
  zero_days      <- matrix( data=1, ncol=max_julian_day, nrow=n_lakeyears )

  # Fill in metadata matrices
  for( ly in 1:n_lakeyears ) {

    # Observations for one particular lake-year
    obs <- filter( imputation_dat, lakeyear==ly )

    # Days on which observations were recorded
    obs_days <- unique( obs$julian_day )

    # Number of days camera data exist for this lake-year
    n_days_per_lakeyear[ly] <- length(obs_days)

    # Days available for imputation
    day_seq              <- 1:n_days_per_lakeyear[ly]
    all_days[ly,day_seq] <- obs_days

    # Observations for one particular lake-year
    imputation_zero_ly <- imputation_zero %>%
                          filter( lakeyear==ly ) %>% 
                          arrange( julian_day )
                          
    imputation_pos_ly  <- imputation_pos %>%
                          filter( lakeyear==ly ) %>%
                          arrange( julian_day )



    # Days of observations
    zero_days_temp <- unique( imputation_zero_ly$julian_day )
    pos_days_temp  <- unique( imputation_pos_ly$julian_day )

    # If there are observations for this particular lake-year, record metadeta
    if( length( zero_days_temp ) > 0 ) {
      zero_lakeyears[ly,zero_days_temp] <- ly
      zero_days[ly,zero_days_temp]      <- zero_days_temp
      zero_nobs[ly,zero_days_temp]      <- imputation_zero_ly %>%
                                           group_by(julian_day) %>% 
                                           summarise(nobs=length(julian_day)) %>%
                                           .$nobs
    }
    if( length( pos_days_temp ) > 0 ) {
      pos_lakeyears[ly,pos_days_temp] <- ly
      pos_days[ly,pos_days_temp]      <- pos_days_temp
      pos_nobs[ly,pos_days_temp]      <- imputation_pos_ly %>%
                                         group_by(julian_day) %>% 
                                         summarise(nobs=length(julian_day)) %>%
                                         .$nobs
    }
  
  }

  if(eval_fit) {
    n_days_per_lakeyear[n_lakeyears+1] <- 0
  }


  #----------------------------------------------------------------------------
  #   Data for covariates
  #----------------------------------------------------------------------------

  false_zero_timeofday      <- false_zero$time_of_day      - 1
  zero_cam_timeofday        <- zero_cam$time_of_day        - 1
  pos_cam_timeofday         <- pos_cam$time_of_day         - 1
  imputation_pos_timeofday  <- imputation_pos$time_of_day  - 1 
  imputation_zero_timeofday <- imputation_zero$time_of_day - 1 

  # List of lake-cameras
  lakecams <- c( false_zero$lake_camera,
                 zero_cam$lake_camera,
                 pos_cam$lake_camera,
                 imputation_zero$lake_camera,
                 imputation_pos$lake_camera )
  lakecams <- sort( unique( lakecams ) )

  # Number of lakecameras
  n_lakecams <- length(lakecams)

  # Update lakecam tags to be between 1 and n_lakecams
  false_zero_lakecams      <- match( false_zero$lake_camera,      lakecams )
  zero_cam_lakecams        <- match( zero_cam$lake_camera,        lakecams )
  pos_cam_lakecams         <- match( pos_cam$lake_camera,         lakecams )
  imputation_pos_lakecams  <- match( imputation_pos$lake_camera,  lakecams )
  imputation_zero_lakecams <- match( imputation_zero$lake_camera, lakecams )

  save( imputation_zero_lakecams, file="imputation_zero_lakecams.Rdata" )

  # Get unique p_seen values
  pseen <- effort_cam %>%
           group_by(lake_camera,p_seen) %>% 
           summarise() %>%
           .$p_seen
  
  #----------------------------------------------------------------------------
  #  Define data and parameters for each model
  #----------------------------------------------------------------------------
  
  if( eval_fit ) {
    ipe <- imputation_pos$ground_count
    ize <- imputation_zero$ground_count
  } else {
    ipe=c(0,0)
    ize=c(0,0)
  }

  jags_data_old <- list( n_lyrs   = n_lakeyears,
                     nlc      = n_lakecams,
                     nzz      = nrow(zero_cam),
                     nzc      = zero_cam$n_zero_cam,
                     nze      = zero_cam$n_zero_eff,
                     nzp      = nrow(false_zero),
                     niz      = nrow(imputation_zero),
                     nip      = nrow(imputation_pos),
                     zp       = false_zero$ground_count,
                     npp      = nrow(pos_cam),
                     ppc      = pos_cam$cam_count,
                     ppe      = pos_cam$ground_count,
                     lc_zz    = zero_cam_lakecams,
                     lc_zp    = false_zero_lakecams,
                     lc_pp    = pos_cam_lakecams,
                     lc_iz    = imputation_zero_lakecams,
                     lc_ip    = imputation_pos_lakecams,
                     tod_zp   = false_zero_timeofday,
                     tod_zz   = zero_cam_timeofday,
                     tod_pp   = pos_cam_timeofday,
                     tod_iz   = imputation_zero_timeofday,
                     tod_ip   = imputation_pos_timeofday,
                     wknd_zp  = false_zero$weekend,
                     wknd_zz  = zero_cam$weekend,
                     wknd_pp  = pos_cam$weekend,
                     wknd_iz  = imputation_zero$weekend,
                     wknd_ip  = imputation_pos$weekend,
                     ly_iz    = imputation_zero$lakeyear,
                     date_iz  = imputation_zero$julian_day,
                     i_iz     = imputation_zero$day_index,
                     ncam_iz  = imputation_zero$n_cam_per_dhly,
                     #ize      = ize,
                     ippc     = imputation_pos$cam_count,
                     ly_ip    = imputation_pos$lakeyear,
                     date_ip  = imputation_pos$julian_day,
                     i_ip     = imputation_pos$day_index,
                     ncam_ip  = imputation_pos$n_cam_per_dhly,
                     #ipe      = ipe,
                     n_days   = n_days_per_lakeyear,
                     all_days = all_days,
                     lyrs_z   = zero_lakeyears,
                     days_z   = zero_days,
                     nobs_z   = zero_nobs,
                     lyrs_p   = pos_lakeyears,
                     days_p   = pos_days,
                     nobs_p   = pos_nobs,
                     pseen_zp = false_zero$p_seen,
                     pseen_zz = zero_cam$p_seen,
                     pseen_pp = pos_cam$p_seen,
                     pseen_iz = imputation_zero$p_seen,
                     pseen_ip = imputation_pos$p_seen,
                     pseen    = pseen )
  
setwd(TMB_dir)
save(jags_data,file="jags_data_old.Rdata")




  # Run all possible models and find model with lowest DIC
  if( run_all ) {

    # Specify locations of model files
    model_file <- list( model_1 = "ImputeEff_h_1.txt",  # no covariates         
                        model_2 = "ImputeEff_h_2.txt",  # cov on alpha            
                        model_3 = "ImputeEff_h_3.txt",  # cov on lambda         
                        model_4 = "ImputeEff_h_4.txt",  # cov on delta            
                        model_5 = "ImputeEff_h_5.txt",  # cov on lambda and alpha
                        model_6 = "ImputeEff_h_6.txt",  # cov on delta and alpha  
                        model_7 = "ImputeEff_h_7.txt",  # cov on delta and lambda
                        model_8 = "ImputeEff_h_8.txt" ) # cov on all             
  
    # Number of models
    n_models <- length( model_file )
  
    # Parameters for each model
    pars <- list( model_1 = c("md","pd","ml","pl","ma","pa","d_track","L","a"),
                               #"MSE_Total","Bias"),
                  model_2 = c("md","pd","ml","pl","ma","pa","mawk","pawk",
                              "d_track","L","a","awk","zdnew","zLnew","zanew",
                              "zawknew"),#"MSE_Total","Bias"),
                  model_3 = c("md","pd","ml","pl","ma","pa","mLwk","pLwk",
                              "d_track","L","a","Lwk","zdnew","zLnew","zLwknew",
                              "zanew"),#"MSE_Total","Bias"),
                  model_4 = c("md","pd","ml","pl","ma","pa","mdtod","pdtod",
                              "d_track","L","a","dtod","MSE_Total","Bias"),
                  model_5 = c("md","pd","ml","pl","ma","pa","mawk","pawk","mLwk",
                              "pLwk","d_track","L","a","awk","Lwk","zdnew",
                              "zLnew","zLwknew","zanew","zawknew","MSE_Total",
                              "Bias"),
                  model_6 = c("md","pd","ml","pl","ma","pa","mdtod","pdtod",
                              "mawk","pawk","d_track","dwd_track","L","a","dtod",
                              "awk","MSE_Total","Bias"),
                  model_7 = c("md","pd","ml","pl","ma","pa","mdtod","pdtod",
                              "mLwk","pLwk","d_track","L","a","dtod","Lwk",
                              "MSE_Total","Bias"),
                  model_8 = c("md","pd","ml","pl","ma","pa","mawk","pawk",
                              "mdtod","pdtod","mLwk","pLwk","d_track","L","a",
                              "awk","dtod","Lwk","MSE_Total","Bias") )
  
    # Initial parameter values for each model
    inits <- list( model_1 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5 )
                             },
                   model_2 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mawk=0, pawk=0.5 )
                             },
                   model_3 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mLwk=0, pLwk=0.5 )
                             },
                   model_4 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mdtod=0, pdtod=0.5 )
                             },
                   model_5 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mawk=0, pawk=0.5, mLwk=0, pLwk=1.5 )
                             },
                   model_6 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mdtod=0, pdtod=0.5, mawk=0, pawk=0.5 )
                             },
                   model_7 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mdtod=0, pdtod=0.5, mLwk=0, pLwk=1.5 )
                             },
                   model_8 = function() {
                               list( md=0, pd=0.5, ml=0, pl=1.5, ma=0, pa=0.5,
                                     mawk=0, pawk=0.5, mdtod=0, pdtod=0.5,
                                     mLwk=0, pLwk=1.5)
                             } )
  
    # Deviance information criterion for each model
    DIC <- numeric( length=n_models )

    # Run each model
    for( i in 1:n_models ) {
  
      cat( "Fitting model", i, "of", n_models, "\n" )
  
      # Run jags
      setwd(input_JAGS_dir)
      jags_output <- jags( data               = jags_data,
                           inits              = inits[[ i ]],
                           parameters.to.save = pars[[ i ]],
                           model.file         = model_file[[ i ]],
                           n.chains           = nchain,
                           n.iter             = nburn+nits,
                           n.burnin           = nburn,
                           n.thin             = nthin,
                           DIC                = TRUE )
  
      # File to save jags output to
      output_file <- paste( "MCMC", i, "txt", sep="." )
  
      # Save jags output
      setwd(output_JAGS_dir)
      save( jags_output, file=output_file )

      # Record DIC
      DIC[i] <- jags_output$BUGSoutput$DIC
  
    }

    # Choose model with lowest DIC to parameterize imputation
    imp_model <- which(DIC==min(DIC))
  
    rm( jags_data, inits, pars, model_file )

  }

  #----------------------------------------------------------------------------
  #  Indices for filling blank camera observations
  #----------------------------------------------------------------------------
  
  # Dataframe for finding similar lakeyears
  lakeyear_data <- imputation_dat %>%
                   distinct( lakeyear, year, area, region, access, launch,
                             camp, resort, dist_kam, dist_van ) %>% 
                   select( lake_ID, lakeyear, year, area, region, access,
                           launch, camp, resort, dist_kam, dist_van ) %>%
                   arrange( lakeyear )

  # Days common to both the index lake and the comparison lakes
  common_days <- array( data = 0,
                        dim  = c( n_lakeyears, n_compl, n_days ) )

  # Days to fill
  uncommon_days <- array( data = 0,
                          dim  = c( n_lakeyears, n_compl, n_days ) )

  # Days for which there are comparison data
  fill_days <- matrix( data = 0,
                       nrow = n_lakeyears,
                       ncol = n_days )

  # Number of days where there's no data for index lake, but there is data for
  # at least one comparison lake
  n_fill <- numeric( length = n_lakeyears )

  # Comparison lakes to use on each particular day
  daily_complake <- array( data = 0,
                           dim  = c( n_lakeyears, n_days, n_compl ) )

  # Number of comparison lakes on each day
  n_daily_complake <- matrix( data = 0,
                              nrow = n_lakeyears,
                              ncol = n_days )

  # Comparison lake-years for each index lake
  comparison_ly <- matrix( data = 0,
                           nrow = n_lakeyears,
                           ncol = n_compl )  
  
  # Number of days in common between lakeyear and comparison lakeyears
  n_common_days <- matrix( data = 0,
                           nrow = n_lakeyears,
                           ncol = n_compl )

  # Number of days where data doesn't exist for lake year but exists for 
  # comparison lakeyears
  n_uncommon_days <- matrix( data = 0,
                             nrow = n_lakeyears,
                             ncol = n_compl )

  # Does comparable lake exist?
  comp_exist <- matrix( data = 0,
                        nrow = n_lakeyears,
                        ncol = n_compl )

  # Number of similar lakeyears for each lakeyear
  n_similar_ly <- numeric( length=n_lakeyears )

  for( ly in 1:n_lakeyears ) {

    similar_lakeyears <- lakeyear_data %>%
                         # Calculate distance of each lake from ly
                         mutate( distance = abs( dist_kam-dist_kam[ly] ) ) %>%
                         # Find similar lakes
                         filter( region   == region[ly],
                                 access   == access[ly],
                                 year     == year[ly],
                                 lakeyear != ly ) %>%
                         # Re-order by distance from lake being filled
                         arrange( distance )

    # Number of similar lakeyears
    n_similar_ly[ly] <- min( nrow(similar_lakeyears), n_compl )

    # If no similar lake exists for this lakeyear, skip to the next lakeyear
    if( n_similar_ly[ly] == 0 )
      next

    # Index sequence for filling in matrices
    ly_seq <- 1:n_similar_ly[ly]

    # Assign most similar lakeyears to comparison lakeyear matrix
    comparison_ly[ly,ly_seq] <- similar_lakeyears[ ly_seq, "lakeyear" ] %>%
                                unlist()

    # Unique julian days in imputation dataset for lakeyear ly
    days_ly <- imputation_dat %>%
               filter( lakeyear==ly ) %>%
               group_by( julian_day ) %>%
               summarise() %>%
               .$julian_day 

    # Loop over similar lakeyears
    for( cly in 1:n_similar_ly[ly] ) {

      # All days in comparison lakeyear
      days_cly <- imputation_dat %>%
                        filter( lakeyear==comparison_ly[ly,cly] ) %>%
                        group_by( julian_day ) %>%
                        summarise() %>%
                        .$julian_day 

      # Days common to lakeyear ly and comparison lakeyear
      common_days_ly <- days_ly[ days_ly %in% days_cly ]

      # Days to be filled by info from comparison lakeyear cly
      uncommon_days_ly <- days_cly[ ! days_cly %in% days_ly ]

      # Number of days common to lakeyear ly and comparison lakeyear
      n_com   <- length( common_days_ly )
      
      # Number of days to be filled by info from comparison lakeyear cly
      n_uncom <- length( uncommon_days_ly )

      if( n_com > 0 && n_uncom > 0 ) {

        # Record that a comparison lake exists
        comp_exist[ ly, cly ] <- 1

        # Number of days common to ly and cly
        n_common_days[ ly, cly ] <- n_com
        
        # Number of days to fill
        n_uncommon_days[ ly, cly ] <- n_uncom

        # Days common to ly and cly
        common_days[ly,cly,1:n_com] <- common_days_ly

        # Days to be filled by in from from cly
        uncommon_days[ly,cly,1:n_uncom] <- uncommon_days_ly

      }

    }  # next comparison-lakeyear

    # If no comparable lake exists for this lakeyear, skip to the next ly
    if( sum( comp_exist[ly, ] == 0 ) )
      next

    # Unique days to be filled in for all comparison-lakeyears
    fill_days_ly <- uncommon_days[ly, , ] %>%
                     as.vector() %>%
                     .[ . != 0 ] %>%
                     unique() %>%
                     sort()

    # Number of unique days to be filled in for all comparison-lakeyears
    n_fill[ly] <- length( fill_days_ly )

    fill_days[ ly, 1:n_fill[ly] ] <- fill_days_ly

    # For each day to be filled, which cly should be used?
    for( t in 1:n_fill[ly] ) {

      day <- fill_days[ly,t]

      cly_to_use <- apply( X      = uncommon_days[ly, , ]==day,
                           MARGIN = 1,
                           FUN    = any)

      n <- length(cly_to_use)

      n_daily_complake[ly,t] <- n

      daily_complake[ly,t,1:n] <- which( cly_to_use )

    }

  } # next lakeyear 

  # Days where there is actual or filled data
  active_days <- matrix( data = 0,
                         nrow = n_lakeyears,
                         ncol = n_days)

  n_active_days <- numeric( length = n_lakeyears )

  for( ly in 1:n_lakeyears ) {

    days_ly <- imputation_dat %>%
                    filter( lakeyear==ly ) %>%
                    group_by( julian_day ) %>%
                    summarise() %>%
                    .$julian_day 

    unique_days <- unique( c( fill_days[ly, ], days_ly ) )

    days <- unique_days[ unique_days>0 ]

    n <- length( days )

    n_active_days[ly] <- n

    active_days[ly,1:n] <- days

  }



  n_comps_per_ly <- rowSums( comp_exist )

  lakeyears_to_fill <- which( n_comps_per_ly > 0 )

  n_lakeyears_to_fill <- length( lakeyears_to_fill )

  comp_lakes <- matrix( data = 0,
                        nrow = n_lakeyears,
                        ncol = n_compl )

  for( ly in lakeyears_to_fill ) {
        
    comp_seq <- 1:n_comps_per_ly[ly]
    
    comp_lakes[ly,comp_seq] <- which( comp_exist[ly, ]==1 )
  
  }

  #----------------------------------------------------------------------------
  #  Impute
  #----------------------------------------------------------------------------

  jags_data <- list( n_lyrs    = n_lakeyears,
                     nlc       = n_lakecams,
                     ncam_iz   = imputation_zero$n_cam_per_dhly,
                     ncam_ip   = imputation_pos$n_cam_per_dhly,
                     nzz       = nrow(zero_cam),
                     nzc       = zero_cam$n_zero_cam,
                     nze       = zero_cam$n_zero_eff,
                     n_days    = n_days_per_lakeyear,
                     nobs_z    = zero_nobs,
                     nobs_p    = pos_nobs,
                     nzp       = nrow(false_zero),
                     zp        = false_zero$ground_count,
                     days_z    = zero_days,
                     days_p    = pos_days,
                     lyrs_z    = zero_lakeyears,
                     lyrs_p    = pos_lakeyears,
                     tod_zp    = false_zero_timeofday,
                     tod_iz    = imputation_zero_timeofday,
                     npp       = nrow(pos_cam),
                     ppc       = pos_cam$cam_count,
                     ppe       = pos_cam$ground_count,
                     wknd_zp   = false_zero$weekend,
                     wknd_zz   = zero_cam$weekend,
                     wknd_pp   = pos_cam$weekend,
                     wknd_iz   = imputation_zero$weekend,
                     wknd_ip   = imputation_pos$weekend,
                     niz       = nrow(imputation_zero),
                     nip       = nrow(imputation_pos),
                     ly_iz     = imputation_zero$lakeyear,
                     ly_ip     = imputation_pos$lakeyear,
                     ippc      = imputation_pos$cam_count,
                     i_iz      = imputation_zero$day_index,
                     i_ip      = imputation_pos$day_index,
                     lc_zz     = zero_cam_lakecams,
                     lc_zp     = false_zero_lakecams,
                     lc_pp     = pos_cam_lakecams,
                     lc_iz     = imputation_zero_lakecams,
                     lc_ip     = imputation_pos_lakecams,
                     date_iz   = imputation_zero$julian_day,
                     date_ip   = imputation_pos$julian_day,
                     all_days  = all_days,
                     nlyfill   = n_lakeyears_to_fill,
                     ncomplks  = n_comps_per_ly,
                     complks   = comp_lakes,
                     comply    = comparison_ly,
                     comp      = common_days,
                     lyfill    = lakeyears_to_fill,
                     ncomp     = n_common_days,
                     t_days    = active_days,
                     nt_days   = n_active_days,
                     npfill    = n_uncommon_days,
                     pfill_day = uncommon_days,
                     nfill     = n_fill,
                     fill_days = fill_days,
                     complyd   = daily_complake,
                     ncomplyd  = n_daily_complake,
                     pseen_zp  = false_zero$p_seen,
                     pseen_zz  = zero_cam$p_seen,
                     pseen_pp  = pos_cam$p_seen,
                     pseen_iz  = imputation_zero$p_seen,
                     pseen_ip  = imputation_pos$p_seen
                     )
  
  parameters <- c("Ely")

  # Load output from model with lowest DIC
  output_file <- paste( "MCMC", imp_model, "txt", sep="." )
  load( output_file )
  attach.jags( jags_output )

  # Use MCMC samples from lowest DIC model to paramaterize initial conditions
  md_mean   <- mean(md)
  pd_mean   <- mean(pd)
  ml_mean   <- mean(ml)
  pl_mean   <- mean(pl)
  ma_mean   <- mean(ma)
  pa_mean   <- mean(pa)
  mLwk_mean <- mean(mLwk)
  pLwk_mean <- mean(pLwk)

  detach.jags()

  inits <- function() {
             list( md   = md_mean,
                   pd   = pd_mean,
                   ml   = ml_mean,
                   pl   = pl_mean,
                   ma   = ma_mean,
                   pa   = pa_mean,
                   mLwk = mLwk_mean,
                   pLwk = pLwk_mean )
           }  

  # Adjust burn-in
  nburn <- 4000


  jags_output <- jags( data               = jags_data,
                       inits              = inits,
                       parameters.to.save = parameters,
                       model.file         = "ImputeEff_h_impute.txt",
                       n.chains           = 2,
                       n.iter             = nburn+nits,
                       n.burnin           = nburn,
                       n.thin             = nthin,
                       DIC                = TRUE )

  save( jags_output, file="Imputation_MCMC.dat" )

  
}
  

posteriors <- function()
{
  load("Model_5_MCMC.dat")
  load("imputation_zero_lakecams.Rdata")
  maxy<-c(25,2,1.3,5,6)
  maxx<-c(1,10,10,5,5)
  minx<-c(0,0,0,0,0)
  panels<-c("A","B","C","D","E")
  xlab=c(expression(paste("P( False-zero camera observation | C"["l,i"]," = 0 ) [logit"^{-1},"(",italic(delta["l"]%*%p["l"])," )]")),
         expression(paste("Mean effort on weekdays | C"["l,i"]," = 0 [exp(",italic(lambda["0,l"]),")]")),
         expression(paste("Mean effort on weekends | C"["l,i"]," = 0 [exp(",italic(lambda["0,l"]),"+",italic(lambda["1,l"]),")]")),
         expression(paste("Effort muliplier on weekdays | C"["l,i"]," > 0 [exp(",italic(alpha["0,l"]),")]")),
         expression(paste("Effort multiplier on weekends | C"["l,i"]," > 0 [exp(",italic(alpha["0,l"]),"+",italic(alpha["1,l"]),")]")))
  attach.jags(jags_output)
  posteriors<-list(d_track,L,Lwk,a,awk)
  new.posts<-list(zdnew,zLnew,zLwknew,zanew,zawknew)
  save(posteriors,file="model posteriors.txt")
  save(new.posts,file="unsampled posteriors.txt")
  
  load("model posteriors.txt")
  load("unsampled posteriors.txt")
  par(mfcol=c(5,1),mar=c(4,1,0.5,1),oma=c(0.5,3,0.5,0.5))
  
  x<-sort(unique(imputation_zero_lakecams))
  x2<-vector()
  for(i in 1:length(x))
  {
    x2[i]<-imputation_zero_lakecams[which(imputation_zero_lakecams==x[i])[1]]
  }
  for(j in 1:5) 
  {
    y<-exp(posteriors[[j]])
    if(j==1)
      y<-1/(1+exp(-posteriors[[j]]))
    if(j==3)
      y<-exp(posteriors[[2]]+posteriors[[3]])
    if(j==5)
      y<-exp(posteriors[[4]]+posteriors[[5]])
    plot(density(y[,1],from=0,to=maxx[j]),ylim=c(0,maxy[j]),main="",xlab=xlab[j],ylab="",cex.lab=1.5,cex.axis=1.2,lwd=1,yaxt="n",col="dark grey")
    for(i in 2:length(posteriors))
    {
      lines(density(y[,i],from=0,to=maxx[j]),lwd=1,col="grey")
    }
    lines(density(new.posts[[j]],from=0,to=maxx[j]),lwd=3)
    legend("topright",panels[j],cex=2,bty="n")
  }
  par(mfcol=c(1,1))
  mtext("Probability density                   ",2,line=2.5,cex=1.5)
}
  
  


  
  
  
  
  
  
  