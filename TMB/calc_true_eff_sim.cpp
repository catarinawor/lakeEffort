#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}



template<class Type>
Type objective_function<Type>::operator()()
{
  
  //likelihood
  DATA_INTEGER(nzz); //number of zero camera observations and zero ground counts - true zeroes
  DATA_INTEGER(nzp); //number of false zero: zero camera effort and positive total effort
  DATA_INTEGER(npp); //number of positive camera observations
  DATA_INTEGER(niz); //number of positive camera observations
  DATA_INTEGER(nip);


  //DATA_INTEGER(zz); //number of zero camera observations and zero ground counts - true zeroes
  //observations for false zeroes -- corresponding ground counts
  //DATA_INTEGER(zp); //ground counts for zero camera observarions
  //DATA_INTEGER(pp); //ground counts for zero camera observarions

  DATA_INTEGER(nlc);

  //observarion for binomial model(true/false zeroes)
  DATA_VECTOR(nzc); //number of zero camera observations (trials)
  DATA_VECTOR(nze); //number of true zeroes (success!)
  
  DATA_VECTOR(zpcount); // number of ground count when 
 

  DATA_VECTOR(ppc); //positive camera counts
  DATA_VECTOR(ppe); // positive ground counts

  //DATA_VECTOR(ncam_iz); //Number of cameras per day-hour-lake-year combinations  
                                            
  
  //indexes or counters
  DATA_IVECTOR(lc_zz); //zero_cam_lakecams -counter for zero camera counts
  DATA_IVECTOR(lc_zp); //false_zero_lakecams - counter for camera counts that were zero
  DATA_IVECTOR(lc_pp); // pos_cam_lakecams  - positive cameracounts

  DATA_IVECTOR(lc_iz);
  DATA_IVECTOR(lc_ip);


  DATA_IVECTOR(ly_iz);
  DATA_IVECTOR(date_iz);
  DATA_IVECTOR(i_iz);
  
  DATA_IVECTOR(ncam_iz);
  DATA_VECTOR(ippc); 
  DATA_VECTOR(ncam_ip);

  DATA_IVECTOR(zero_agg_ind); 
  DATA_INTEGER(lzero_agg_ind);
  
  DATA_VECTOR(pos_agg_ind);
  DATA_INTEGER(lpos_agg_ind);

  DATA_VECTOR(allEff_ind_zero); 
  DATA_VECTOR(allEff_ind_pos); 
  DATA_INTEGER(nallEff_ind); 

  

  
  //PARAMETER(d_track);
  //hyper parameter
  PARAMETER(md);
  PARAMETER(pd);

  PARAMETER(ml);
  PARAMETER(pl);

  PARAMETER(ma);
  PARAMETER(pa);

  //need to define these
  
  PARAMETER_VECTOR(d); //logit transformed probabilities of false zero
  PARAMETER_VECTOR(L); // mean value for false zero obs in log space
  PARAMETER_VECTOR(a); // multiplier of camera counts to match ground counst

  Type ans=0.;

  Type sdd;
  Type sdl;
  Type sda;


   //probability of true zeroes intermediate values
  //vector<Type> prob(zz);
  //vector<Type> delta(zz);

  vector<Type> prob(nzz);
  vector<Type> delta(nzz);

  //false negatives intermediate values
  
   vector<Type> lambda(nzp);

  //true positives intermediate values
  //vector<Type> alpha(pp);
  vector<Type> alpha(npp);
  vector<Type> lambpp(npp);


  //imputation intermediate values
  vector<Type> pos(niz);
  vector<Type> delta_iz(niz);
  vector<Type> prob_iz(niz);
  vector<Type> lambda_iz(niz);
  vector<Type> eff(niz);
  vector<Type> Iz(niz);
  vector<Type> aggIz(lzero_agg_ind);
  aggIz.setZero();

  vector<Type> meff(nip); 
  vector<Type> alpha_ip(nip);
  vector<Type> Iptemp(nip);
  vector<Type> Ipfull(nip);
  vector<Type> Ip(nip);
  vector<Type> aggIp(lpos_agg_ind);
  aggIp.setZero();
  
  vector<Type> sumEff(nallEff_ind);
  sumEff.setZero();


  //Type ffpen=Type(0.);
  //Type dummy =Type(0);
    
  
  //hyper priors for the probability of a true zero effort observation 
  // uniform prior on mean (logit transformed) probability of true zeroes
  //dummy = posfun(Type(md),Type(-5),ffpen);   //  added to keep parameter within bounds
  //Type tmpp = Type(5)-Type(md);
  //dummy = posfun(tmpp,Type(0.0),ffpen);   // added to keep parameter within bounds       
  //ans += -log(Type(1.)/(Type(5)-Type(-5)));
  ans -= dnorm(md,Type(0.0), Type(2.0),1);
  //ans += -ffpen;
  // gamma prior on precision probablitity of true zeroes
  ans -= dgamma(pd,Type(0.001), Type(1000),1);       // Assign N(0,1) distribution u 
  sdd = sqrt(Type(1.0)/pd);

  // hyper prior for count of a positive effort observation given zero camera observation
  //uniform prior on mean of true value when zero obs are seen by camera
  //dummy= posfun(Type(ml),Type(-5),ffpen);   // BvP added to keep parameter within bounds
  //Type tmp = Type(5)- Type(ml);
  //dummy = posfun(tmp,Type(0.0),ffpen);   // BvP added to keep parameter within bounds        
  //ans += -log(Type(1.)/(Type(5)-Type(-5)));
  ans -= dnorm(ml,Type(0.0), Type(2.0),1);
  //ans += -ffpen;
  // gamma prior on precision
  ans -= dgamma(pl,Type(0.001), Type(1000),1);       // Assign N(0,1) distribution u 
  sdl = sqrt(Type(1.0)/pl);

  //hyper priors for the magnitude of a positive effort observation given positive camera observation  
  // uniform prior on log mean distribution of positive camera counts
  //dummy = posfun(Type(ma),Type(-5),ffpen);   // added to keep parameter within bounds
  //Type tmppp = Type(5)-Type(ma);
  //dummy = posfun(tmppp,Type(0.0),ffpen);   //  added to keep parameter within bounds       
  ans -= dnorm(ma,Type(0.0), Type(2.0),1);
  //ans += -log(Type(1.)/(Type(5)-Type(-5)));
  //ans += -ffpen;
  ans -= dgamma(pa,Type(0.001), Type(1000),1);       // Assign N(0,1) distribution u 
  sda = sqrt(Type(1.0)/pa);

  //priors
  //# includes lake-cameras that need imputation but have no recorded effort to go with camera effort
  for(int i=0;i<nlc;i++){

      //Prior in d 
      ans += -dnorm(d(i),md,sdd,1);
    
      //Prior in L
      ans += -dnorm(L(i),ml,sdl,1); 

      //Prior in a
      ans += -dnorm(a(i),ma,sda,1);
  

  }

 
    //probability of true zeroes
    for(int y=0;y<nzz;y++){

     //Prior in d 
     //ans += -dnorm(d(lc_zz(y)),md,sqrt(Type(1.0)/pd),true);
   
      //likelihood of d   
      //delta(y) =Type(1.0)/(Type(1.0)+exp(-d(lc_zz(y))));  //logit back transform delta
      delta(y) = invlogit(d(lc_zz(y)));
      prob(y)= Type(1.0)-delta(y);//probability of seeing a true zero
      //ans+= - dbinom(nze(y),nzc(y),prob(y),1); //binomial likelihood
      ans+= - dbinom(nze(y),nzc(y),prob(y),1); //binomial likelihood
      

    }

    //counts when camera obs is zero - false negatives
    for(int b=0;b<nzp;b++){

      //Prior in L
      //ans += -dnorm(L(lc_zp(b)),ml,sqrt(Type(1.0)/pl),true); 

      //Likelihood
      lambda(b) = exp(L(lc_zp(b))); //mean value for false zero obs 
      ans+= -dpois(zpcount(b),lambda(b),1);  // poisson likelihood for false zeroes (actual numbes that should have been seen)
      
    }

    //true positives
    for(int i=0;i<npp;i++){

      //Prior in a
      //ans += -dnorm(a(lc_pp(i)),ma,sqrt(Type(1.0)/pa),true);
  
      alpha(i)= exp(a(lc_pp(i)));
      lambpp(i)= ppc(i)*alpha(i); // mean distribution of positive camera observations
      ans+= -dpois(ppe(i),lambpp(i),1);

    }

    SIMULATE {

    // imputations for zero camera effort
    for(int ii=0;ii<niz;ii++){
        
      delta_iz(ii) = invlogit(d(lc_iz(ii)));//*pseen_iz[i]
    
      prob_iz(ii)= Type(1.0)-delta_iz(ii);
    
      pos(ii) = rbinom(Type(1.0),prob_iz(ii));
      
      lambda_iz(ii) = exp(L(lc_iz(ii)));
      
      eff(ii) = rpois(lambda_iz(ii));  // use lambda to predict a positive effort given zero camera effort
    
      Iz(ii) = pos(ii) * eff(ii)/ncam_iz(ii);
      
      aggIz(zero_agg_ind(ii)) += Iz(ii);
      
    }

    for(int bb=0;bb<nip;bb++){
      
      alpha_ip(bb) = exp(a(lc_ip(bb)));
     
      meff(bb) = alpha_ip(bb)*ippc(bb);   // predict mean effort using alpha and observed camera effort
     
      Iptemp(bb) = rpois(meff(bb));
     
      Ipfull(bb) = Iptemp(bb)/ncam_ip(bb);
      
      //Ip(bb) = Ipfull(bb);      // predict effort for a positive camera observation
      
      //aggIp(pos_agg_ind(bb)) += Ipfull(bb);

    
    }

    //sum effort for each lakeyear and for each julian day

    //for(int n=0;n<nallEff_ind;n++){

     // std::cout<<"aqui"<<std::endl;

        //if(isNA(allEff_ind_pos(n)) ){
        //   sumEff(n) += aggIp(allEff_ind_pos(n));
        //} else{
        //  sumEff(n) += aggIp(allEff_ind_pos(n));
        //}

        //if( !isNA(allEff_ind_zero(n)) ){
        // sumEff(n) += aggIz(allEff_ind_zero(n));
        //}}

    
    
    

      REPORT( eff);
      REPORT( Iz);
      REPORT( Ipfull);
      REPORT( aggIp);
      REPORT( aggIz);
      REPORT( sumEff);

    }

   // Imputation of effort  --------------------------------------------------------------------------------------------------------------
   //    for(i in 1:niz)   # imputations for zero camera effort
   //    {
   //       pos[i]~dbern(delta_iz[i])   # use delta to predict a one or a zero given zero camera effort
   //       logit(delta_iz[i])<-d[lc_iz[i]]*pseen_iz[i]
   //       eff[i]~dpois(lambda_iz[i])  # use lambda to predict a positive effort given zero camera effort
   //       log(lambda_iz[i])<-L[lc_iz[i]]
   //       Iz[ly_iz[i],date_iz[i],i_iz[i]]<-pos[i]*eff[i]/ncam_iz[i]   # delta poisson distribution
   

   //        MSE1[i]<-pow(Iz[ly_iz[i],date_iz[i],i_iz[i]]-ize[i],2.)
   //    }



  //  }

    //imputation of values

    // imputations for zero camera effort
    //for(int yy=0;yy<niz;yy++){

    //  delta_iz(yy) = (Type(1.0)/(Type(1.0)+exp(-d(lc_iz(yy))))) *pseen_iz(yy);
    //  pos(yy) = rbinom(delta_iz(yy));

    //  lambda_iz(yy) = exp(L(lc_iz(yy)));
    //  eff((yy)) = rpois(lambda_iz(i))  // use lambda to predict a positive effort given zero camera effort

    //  Iz[yy] = pos(yy)*eff(yy)/ncam_iz(yy)

    //}


    //for(i in 1:niz)   # imputations for zero camera effort
    //   {
    //      pos[i]~dbern(delta_iz[i])   # use delta to predict a one or a zero given zero camera effort
    //      logit(delta_iz[i])<-d[lc_iz[i]]*pseen_iz[i]
    //      eff[i]~dpois(lambda_iz[i])  # use lambda to predict a positive effort given zero camera effort
    //      log(lambda_iz[i])<-L[lc_iz[i]]
    //      Iz[ly_iz[i],date_iz[i],i_iz[i]]<-pos[i]*eff[i]/ncam_iz[i]   # delta poisson distribution
    //      
    //   }
   
    //   for(i in 1:nip)   # imputations for positive camera effort
    //   {
    //      log(alpha_ip[i])<-a[lc_ip[i]]
    //      meff[i]<-alpha_ip[i]*ippc[i]   # predict mean effort using alpha and observed camera effort
    //      Iptemp[i]~dpois(meff[i])   # predict effort for a positive camera observation
    //      Ipfull[i]<-Iptemp[i]/ncam_ip[i]
    //      #logit(alpha_ip[i])<-a[lc_ip[i]]
    //      #Iptemp[i]~dnegbin(alpha_ip[i],ippc[i])
    //      #Ipfull[i]<-(Iptemp[i]+ippc[i])/ncam_ip[i]
    //      Ip[ly_ip[i],date_ip[i],i_ip[i]]<-Ipfull[i]       # predict effort for a positive camera observation
    //      #     MSE2[i]<-pow(Ip[ly_ip[i],date_ip[i],i_ip[i]]-ipe[i],2.)
    //   }

    //   for(i in 1:n_lyrs)
    //  {
    //      Iz[i,1,1]<-0
    //      Ip[i,1,1]<-0
    //      for(k in 1:n_days[i])
    //      {
    //         sumE[i,all_days[i,k]]<-sum(Iz[lyrs_z[i,all_days[i,k]],days_z[i,all_days[i,k]],1:nobs_z[i,all_days[i,k]]])+
    //                       sum(Ip[lyrs_p[i,all_days[i,k]],days_p[i,all_days[i,k]],1:nobs_p[i,all_days[i,k]]])
    //      }
    //   }
       


    //ADREPORT(lambpp);
    //ADREPORT(lambda);
    //REPORT(delta);
    REPORT(md);
    REPORT(pd);
    REPORT(ma);
    REPORT(pa);
    REPORT(ml);
    REPORT(pl);
    //REPORT(prob);
    //REPORT(nzz);
    //REPORT(lambda);
    REPORT(d);
    REPORT(L);
    //REPORT(ffpen);
    REPORT(a);
    //REPORT(alpha);
    //ADREPORT(prob);
    
    return ans; 


}



