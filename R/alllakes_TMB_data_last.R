################################################################################
# imputeEffort.R
# Author: Catarina Wor, based on code by
# Steven Rossi (spr1@sfu.ca)
# June, 2018
#
#
# Input files: (a) "Camera Effort2.data"
#              (b) "Lake Info.data"
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
################################################################################
# Load necessary packages
library(dplyr)

#this line must be changed if ran in another machine
source("~/Documents/lake_effort/R/directories.R")
setwd(input_data_dir)

#----------------------------------------------------------------------------
#  Build effort and lake datasets
#----------------------------------------------------------------------------
# Build camera effort dataset
effort_dato <- read.table( "Camera Effort2.data", header=TRUE ) %>%
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
  
  #summary(effort_dato )

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
  #summary(lake_dat)

  # Get the individual lake names from the effort data set
  lake_names <- as.character(unique(effort_dato$name))
 
  # Identify lakes in effort dataset by numeric ID codes
  effort_dato$lake_ID <-as.numeric(effort_dato$name)
  
  # Identify lake information for each observation in the effort dataset
  lake_dat_obs <- lake_dat[ effort_dato$lake_ID, ]
  
  # Merge the effort dataset with the corresponding lake information, then
  # delete unwanted/invalid observations
  effort_dat <- cbind(effort_dato,lake_dat_obs) %>%
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


    rm( lake_names,
      lake_dat,
      lake_dat_obs )
  #summary(effort_dat)
  #----------------------------------------------------------------------------
  #  Partition effort dataset by method of observation
  #----------------------------------------------------------------------------

  # Ground observations
  effort_ground <- effort_dat %>%
                   filter( method=="GR" ) %>%
                   rename( ground_count = angler_count ) %>%
                   select( -method )
  summary(effort_ground)
  # Camera observations
  effort_cam <- effort_dat %>%
                filter( method=="CAM" ) %>%
                mutate( lake_camera = 10*lake_ID+location) %>%
                rename( cam_count = angler_count ) %>%
                select( -method )
  summary(effort_cam)
  
  sort(unique(effort_ground$lake_ID))
  sort(unique(effort_cam$lake_ID))


  #----------------------------------------------------------------------------
  #  Create new effort dataset with ground and cam observations on same line
  #----------------------------------------------------------------------------

  # Create new dataset to hold both camera and ground observations
  # Start off with the same column names as the original effort dataset
  effort <- effort_dat[0, ]

  # Get unique locations from effort dataset
  locations <- unique(effort_dat$location)

  # Loop over all locations. At each location, combine the camera observations
  # with their corresponding ground observations.
  for( i in locations ) {

    # Camera observations at location i
    effort_cam_i <- filter( effort_cam, location==i )
    
    #length(unique(effort_ground$datetimelake))
    #length(unique(effort_cam_i$datetimelake))

    # Find the ground observations that correspond to camera observations i
    cam_index <- match( effort_ground$datetimelake, effort_cam_i$datetimelake )
    match( effort_ground$datetimelake, effort_cam$datetimelake )

    length( cam_index)
    # Combine ground and camera observations
    effort_i <- cbind( effort_ground,
                       cam_count   = effort_cam_i$cam_count[cam_index],
                       lake_camera = effort_cam_i$lake_camera[cam_index] )

    # Add observations at location i to effort dataset
    effort <- rbind( effort, effort_i )
    
  }

  
  unique(effort_ground$lakeyear)
  unique(effort_cam$lakeyear)
  unique(effort$lakeyear)

  #Effort only contain observations that are present for both ground and camera

  # Convert effort matrix to dataframe and delete datetimelake column
  effort <- effort %>%
            as.data.frame() %>%
            select( -datetimelake )

  
  #----------------------------------------------------------------------------
  #  Partition effort dataset by type of observations (Positive, zero, etc...)
  #----------------------------------------------------------------------------

  # Positive camera observations -- when ground counts were present
  pos_cam <- filter( effort, cam_count>0 )
  summary(pos_cam)


  # Zero camera observations and zero ground counts -- when ground counts were present
  # n_zero_cam -> number of zero obs by camera in lake,date and time
  # n_zero_eff -> number of true zeroes -ground obs in lake,date and time
  zero_cam <- effort %>%
              group_by( area, weekend, time_of_day, hour, year, #p_seen
                        lake_camera, lake_ID,lakeyear ) %>%
              arrange(  area, weekend, time_of_day, hour, year, #p_seen
                        lake_camera, lake_ID,lakeyear  ) %>%
              filter( cam_count==0 ) %>%
              mutate( ground_count=ifelse( test = ground_count>0,
                                           yes  = 1,
                                           no   = 0 ) ) %>%
              summarise( n_zero_cam = sum(cam_count==0),
                         n_zero_eff = sum(ground_count==0) ) 
  
  head(as.data.frame(zero_cam))
  #summary(zero_cam)
  #dim(zero_cam)
  # False zero: zero camera effort and positive total effort
  false_zero <- effort %>%
                filter( ground_count > 0,
                        cam_count   == 0 )
  
                

  lakeyears <- sort(unique(effort_cam$lakeyear))

  # Number of days per lakeyear
 


  dim(false_zero)
  summary(false_zero)

  lakecams <- c( false_zero$lake_camera,
                 zero_cam$lake_camera,
                 pos_cam$lake_camera )

  lakecams_u <- sort( unique( lakecams ) )
  lakecams_pu<-sort(unique(pos_cam$lake_camera))

  # Number of lakecameras
  n_lakecams <- length(lakecams_u)

  

  # Update lakecam tags to be between 1 and n_lakecams
  false_zero_lakecams      <- match( false_zero$lake_camera,      lakecams_u )
  zero_cam_lakecams        <- match( zero_cam$lake_camera,        lakecams_u )
  pos_cam_lakecams         <- match( pos_cam$lake_camera,         lakecams_pu )

(match( pos_cam$lake_camera,         lakecams_u ))

pos_cam  #<- Positive camera observations -- when ground counts were present


summary(false_zero)
summary(zero_cam) #n_zero_cam = sum(cam_count==0)
                  #n_zero_eff = sum(ground_count==0)
summary(pos_cam)

sort(unique(false_zero$lake_camera))
sort(unique(zero_cam$lake_camera))
order(sort(unique(pos_cam$lake_camera)))

#the imputation dat has more lakecams because the parameters are random
sort(unique(false_zero_lakecams -1  )   )
sort(unique(zero_cam_lakecams     -1 )  )
sort(unique(pos_cam_lakecams       -1 ) )

# Get unique p_seen values
#  pseen <- effort_cam %>%
#           group_by(lake_camera,p_seen) %>% 
#           summarise() %>%
#           .$p_seen
#
#----------------------------------------------------------------------------
#  TMB data
#----------------------------------------------------------------------------
  

  TMB_data<- list(
                #likelihood counters
                  nzz      = nrow(zero_cam),
                  nzp      = nrow(false_zero),
                  npp      = nrow(pos_cam),

                  zz      = length(unique(zero_cam_lakecams)),
                  zp      = length(unique(false_zero_lakecams)),
                  pp      = length(unique(pos_cam_lakecams)),

                  nzc      = zero_cam$n_zero_cam, #sum(cam_count==0)
                  nze      = zero_cam$n_zero_eff, #sum(ground_count==0)
                  zpcount  = false_zero$ground_count, #real ground counts when camera saw no effort
                  ppc      = pos_cam$cam_count, #real observed effort
                  ppe      = pos_cam$ground_count, #effort seen by camera
                  lc_zz    = zero_cam_lakecams-1, #counter for zero camera counts
                  lc_zp    = false_zero_lakecams-1,
                  lc_pp    = pos_cam_lakecams-1
                  

                  )



setwd(TMB_dir)
save(TMB_data,file="TMB_data.Rdata")









  #----------------------------------------------------------------------------
  #  Imputation data
  #----------------------------------------------------------------------------
   cutoff    = 0.3 #cuttof to impute data, cameta has to bserve at least 30% of lake area

  # Imputation dataset: Positive camera effort but no ground counts
  imputation_dat <- effort_cam %>%
                    filter( year >= 2009 ) %>%
                    mutate( ground_count = NA ) %>%
                    arrange( lake_ID )
  summary(  imputation_dat )

   # Total number of days in imputation dataset
  n_days_total <- length(unique(imputation_dat$julian_day))


  # Unique lakeyears
  lakeyears <- effort_cam %>%
               group_by( lakeyear ) %>%
               summarise() %>%
               .$lakeyear
  


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

  lakeyears <- sort(unique(imputation_dat$lakeyear))

  # Number of lakeyears
  n_lakeyears <- length( lakeyears )

#----------------------------------------------------------------------------
#  Partition imputation dataset and record metadata
#----------------------------------------------------------------------------


  imputation_dat <- arrange( imputation_dat, lakeyear )

  # Update lakeyear tags to range from 1 to n_lakeyears
  imputation_dat$lakeyear <- match(imputation_dat$lakeyear,lakeyears)
  false_zero$lakeyear     <- match(false_zero$lakeyear,lakeyears)
  pos_cam$lakeyear        <- match(pos_cam$lakeyear,lakeyears)


  summary(imputation_dat$lakeyear)
  
  # Create codes to differentiate observations
  imputation_dat <- imputation_dat %>%
                    mutate( lakeyearday     = 1000*lakeyear +
                                                   julian_day,
                            dayhourlakeyear = 100000*julian_day +
                                                1000*hour +
                                                     lakeyear )



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
  all_days <- matrix( data=1, nrow=n_lakeyears, ncol=n_days)

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
                 pos_cam$lake_camera )#,
                 #imputation_zero$lake_camera,
                 #imputation_pos$lake_camera )


  lakecams_u <- sort( unique( lakecams ) )

  # Number of lakecameras
  n_lakecams <- length(lakecams_u)

  # Update lakecam tags to be between 1 and n_lakecams
  false_zero_lakecams      <- match( false_zero$lake_camera,      lakecams_u )
  zero_cam_lakecams        <- match( zero_cam$lake_camera,        lakecams_u )
  pos_cam_lakecams         <- match( pos_cam$lake_camera,         lakecams_u )
  imputation_pos_lakecams  <- match( imputation_pos$lake_camera,  lakecams_u )
  imputation_zero_lakecams <- match( imputation_zero$lake_camera, lakecams_u )


imputation_zero_filter<-subset(imputation_zero,  imputation_zero_lakecams>0)
imputation_pos_filter<-subset(imputation_pos,  imputation_pos_lakecams>0)

imputation_pos_lakecams  <- match( imputation_pos_filter$lake_camera,  lakecams_u )
imputation_zero_lakecams <- match( imputation_zero_filter$lake_camera, lakecams_u )


sort(unique(false_zero$lake_camera))
sort(unique(zero_cam$lake_camera))
sort(unique( pos_cam$lake_camera))
sort(unique(imputation_pos$lake_camera))
sort(unique(imputation_zero$lake_camera))



summary(aa)
sort(unique(aa$lake_camera))


#the imputation dat has more lakecams because the parameters are random
sort(unique(false_zero_lakecams -1  )   )
sort(unique(zero_cam_lakecams     -1 )  )
sort(unique(pos_cam_lakecams       -1 ) )
sort(unique(imputation_pos_lakecams -1 ))
sort(unique(imputation_zero_lakecams -1))



# Get unique p_seen values
  pseen <- effort_cam %>%
           group_by(lake_camera,p_seen) %>% 
           summarise() %>%
           .$p_seen
  


#----------------------------------------------------------------------------
#  TMB data
#----------------------------------------------------------------------------
  


  TMB_data<- list(
                #likelihood counters
                  nzz      = nrow(zero_cam),
                  nzp      = nrow(false_zero),
                  npp      = nrow(pos_cam),
                #imputation counters
                  niz      = nrow(imputation_zero)
            
            #integer vectors
                  nzc      = zero_cam$n_zero_cam,
                  nze      = zero_cam$n_zero_eff,
                  zp       = false_zero$ground_count,
                  ppc      = pos_cam$cam_count,
                  ppe      = pos_cam$ground_count,
                  lc_zz    = zero_cam_lakecams-1,
                  lc_zp    = false_zero_lakecams-1,
                  lc_pp    = pos_cam_lakecams-1
                  

                  )


setwd(TMB_dir)
save(TMB_data,file="TMB_data.Rdata")

  
  
  
  
  