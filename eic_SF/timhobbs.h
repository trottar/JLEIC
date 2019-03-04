/* This is the collection of subroutine from Tim. Hobbs's pion-nucleon
 * fracture fuction to obtain the pion structure function
 * 
 * The original code is writtne by FORTRAN 
 * Cpp conerting by K.Park
 * 
 * Motivation:
 * Changing the pion FF with various assumption based on the pion exchange model
 * J. McKenney, N. Sato, W. Melnitchouk PRD 93, 054011 (2016) 
*/





//
//   Exactly quoted from Tim's Fortran code.
//
//WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
//      L = 1.18D0      !COV DIPOLE: HSS NORM
//      L = 1.71D0      !IMF DIPOLE: HSS NORM
//      L = 1.63D0       !HSS +
const double      L = 1.560;      //!HSS CENT. VALUE
//      L = 1.48D0      !HSS -
//C***********************************************************************
//!HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
//!      typ = 1 !DIPOLE FORM FACTOR
const int typ = 2; // !EXPONENTIAL FORM FACTOR
//!      typ = 3 !COV. DIPOLE FORM FACTOR
//C***********************************************************************
//!THE PARAMETER 'dis' SELECTS THE CHARGE CONFIG. OF THE INTERM. STATE
//!       dis = 0 !CHARGE EXCHANGE
const int dis = 1; // !NEUTRAL EXCHANGE
//C***********************************************************************
//!THIS CHOOSES AMONG THE VARIOUS POSSIBLE COMBINATIONS OF SPLITTING FUNCTION/PDF
const int FLAG = 0;
//!THE FLAGS TOGGLE AMONG DISSOCIATION MODES AS: 
//!     FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
//!     FLAG = 1  --- THE RHO CONTRIBUTION       | J = 1 + 1/2
//     FLAG = 2  --- THE PI-DELTA CONTRIBUTION  | J = 0 + 3/2
//!---------------------------------------------------------------------




// ******************************************************************************
/* FUNCTION GIVING f(y) FOR N-PION VERTEX, WHERE
 *    y IS THE IMF MOMENTUM OF THE INTERMEDIATE MESON.
 *
 *  WRITTEN: W. Melnitchouk (1999)
 *  MODIFIED: T. HOBBS (2013)
 */
// ******************************************************************************

double fypiN(double y,double kT,double L,int typ,int dis){

  double  ss,kT2,SpiN;
  double fypiN,t;
  double pi,mN,mpi=0,mP=0,g_piNN,gg=0,FF=0,sM;

  pi = 4*atan(1.0);
  mN  = 0.93891897;  //!masses in GeV!!
  if(dis==0){
    //!!**** THE DISSOCIATION P --> Dbar0 + LAMBDA_c^+ 
    mpi = 0.13957018;   // !THE MASS OF THE pi^-
    mP = mN;            //  !THE MASS OF THE PROTON
    //!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
    g_piNN = sqrt(14.40 * 4*pi);                 //! g_{pi NN}
    gg = 2.0 * pow(g_piNN,2.) / (16.0 * pi*pi);   //!2 ISOSPIN FACTOR
  }
  else if(dis==1){
    //!!**** THE DISSOCIATION P --> Dbar0 + SIGMA_c^+ 
    mpi = 1.865;   //!THE MASS OF THE Dbar0
    mP = 2.4529;   //!THE MASS OF THE CHARMED SIGMA
    //!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
    g_piNN = sqrt(0.5760 * 4*pi);           //! as N-D-Sigma_c
    gg = pow(g_piNN,2.) / (16.0 * pi*pi);   //!1 ISOSPIN FACTOR
  }
  if(y<=0.0 || y>=1.0){
    fypiN = 0.;
    return fypiN;
  }
  kT2 = kT*kT;
  SpiN = (kT2 + mpi*mpi)/y + (kT2 + mP*mP)/(1.0-y);

  if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SpiN));             //! monopole
  else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + SpiN)),2.) ;     //! dipole
  else if(typ==2)  FF = exp( (mN*mN - SpiN)/(L*L) );                //! expon
  else if(typ==3) {
    t = (- kT2 - mN*mN*y*y) /(1.0-y);
    FF = pow(((L*L - mpi*mpi) / (L*L - t)),2.);                      //! cov dip
  }
  else if(typ==4){                       
    sM = (kT2 + (1.0+y)*mpi*mpi)/y +
      (kT2 + y*mP*mP)/(1.0-y) + mN*mN;
    FF = (pow(L,4.) + pow(mP,4.))/(pow(L,4.) + sM*sM); // ! DIPOLE -- s-channel Lambda exchange
  }
  
  ss = ( kT2 + pow((mP - (1.0-y) * mN),2.)) / (1.0-y)
    / pow(( (1.0-y)*(SpiN - mN*mN) ),2.) * FF*FF;

  fypiN =  gg * (1.0-y) / y * ss ;
  return fypiN;

}














// ******************************************************************************
/*  Function giving numerical value of f(y) for N-rho-N
 *  OUTPUT IS THE NEUTRON ---> rho^- + PROTON SPLITTING FUNCTION
 *  By: T. Hobbs on NOV 12, 2013
 *  taken from notes "spin-1, m_B /= m_N," May 17, 2012
 ****************************************************************************
 */
double f_rhoN(double y,double kT,double L,int typ,int dis){

  double kT2,SRoN,P_k,pl_k,P_p,t,sv,st,si,ss;
  double f_rhoN,sM;
  double pi,mN,mrho,mP=0,g_RoNN,f_RoNN,gg=0,fff=0,fg=0,FF=0;

  pi = 4*atan(1.0);
  mN  = 0.93891897; // !masses in GeV!!

  if(dis==0){
    //!**** THE DISSOCIATION N --> rho^- + PROTON
    mrho = 0.7754;   //!THE MASS OF THE rho^-
    mP = mN;         //!THE MASS OF THE PROTON
    //!*** WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
    g_RoNN = sqrt (2.0 * 0.55 * 4. * pi);   //!rho-N-N (Hohler and Pieteranin)
    f_RoNN = 6.10 * g_RoNN;
    gg = pow(g_RoNN,2.) / (16.0 * pi*pi);
    fff = pow(f_RoNN,2.) / (16.0 * pi*pi);
    fg = f_RoNN*g_RoNN / (16.0 * pi*pi);    // !FACTOR OF 2: ISOSPIN
  }
  else if(dis==1){
    //!**** THE DISSOCIATION P --> Dbar*0 + LAMBDA_c^+ 
    mrho = 0.77540;     //!THE MASS OF THE rho^-
    mP = mN;            //!THE MASS OF THE PROTON
    //!*** WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
    g_RoNN = sqrt (0.550 * 4.0 * pi);    //  !rho-N-N (Hohler and Pieteranin)
    f_RoNN = 6.10 * g_RoNN;
    gg = pow(g_RoNN,2.) / (16.0 * pi*pi);
    fff = pow(f_RoNN,2.) / (16.0 * pi*pi);
    fg = f_RoNN*g_RoNN / (16.0 * pi*pi);
  }

  if(y<=0.0 || y>=1.0){
    f_rhoN = 0.0;
    return f_rhoN;
  }

  kT2 = kT*kT;
  SRoN = (kT2 + mrho*mrho)/y + (kT2 + mP*mP)/(1.0-y);
     
  if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SRoN));           //! MONOPOLE
  else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + SRoN)),2.);    //! DIPOLE
  else if(typ==2)  FF = exp( (mN*mN - SRoN)/(L*L) );              //! EXPONENTIAL
  else if(typ==3) {
    t = (- kT2 - mN*mN*y*y) /(1.0-y);
    FF = pow(((L*L - mrho*mrho) / (L*L - t)),2.);                  //! COV. DIPOLE
  }
  else if(typ==4){ 
            sM = (kT2 + (1.0+y)*mrho*mrho)/y +
	      (kT2 + y*mP*mP)/(1.0-y) + mN*mN;  // ! DIPOLE -- s-CHANNEL LAMBDA EXCHANGE
            FF = (pow(L,4.) + pow(mP,4.))/(pow(L,4.) + sM*sM);
  }
         
  //C...TOPT WITH P[alpha] - p[alpha] DERIVATIVE COUPLING
  P_k = (mrho*mrho + y*y*mN*mN + kT2)/2.0/y;
  P_p = (mP*mP + pow((1.0-y),2.)*mN*mN + kT2)/2.0/(1.0-y);
  pl_k = (mP*mP+kT2)*y/2.0/(1.0-y) 
    + (mrho*mrho+kT2)*(1.0-y)/2.0/y + kT2;
  
  sv = -6.0*mN*mP + 4.0*P_k*pl_k/(mrho*mrho) + 2.0*P_p;
  st = -(pow(P_p,2.)) + P_p*pow((mP+mN),2.) - mP*mN*(mP*mP+mN*mN+mP*mN)
    + 1.0/(2.0*mrho*mrho) * ( (P_p - mP*mN)*pow((P_k-pl_k),2.)
			      - 2.0*(P_k-pl_k)*(mP*mP*P_k - mN*mN*pl_k)
			      + 2.0*P_k*pl_k*(2.0*P_p-mP*mP-mN*mN) ); 

  si = -4.0*(mP+mN)*(mP*mN - P_p) 
    -2.0*(mP*P_k*P_k - (mP+mN)*P_k*pl_k + mN*pl_k*pl_k)/(mrho*mrho);
          
  if(kT>0.0 && (sv<0.0 || st<0.0)){
    cout << "CS1 -- ##### kT,y,sv,st =" << kT << ", " << y <<", "<< sv <<", "<< st<< endl;	    
	   // should be stopped !!
  }
        
  ss = (gg*sv + fff*st/mN*mN + fg*si/mN)    //!THE FULL
    / pow(( y*(SRoN - mN*mN) ),2.) * FF*FF;      //!EXPRESSION (UN-INT.)**

  f_rhoN = y / (1.0-y) * ss;

  return f_rhoN;

}











// ******************************************************************************
/*
 *  Function giving numerical value of f(y) for N-Rho-Delta
 *  OUTPUT IS THE NEUTRON ---> RHO + DELTA SPLITTING FUNCTION
 *  By: T. Hobbs on JAN 7, 2014
*/
// ******************************************************************************

double f_RhoDel(double y,double kT,double L,int typ,int dis){

  double  kT2,SRoD,P_k,pl_k,P_p,t;
  double f_RhoDel,sr,ss;
  double pi,mN,mD,mRo,g_NDS,gg=0,FF=0,sM;
     
  pi = 4*atan(1.);
  mN  = 0.93891897; //!masses in GeV!!
    if (dis==0){
      //!!**** THE DISSOCIATION N --> (RHO^0 + DELTA^0) + (RHO^- + DELTA^+)
      mRo = 0.7754;    //!THE MASS OF THE RHO MESON
      mD = 1.232;      // !Mass OF THE DELTA ISOBAR
      //!!*** WE USE THE COUPLINGS INFERRED FROM HOLZENKAMP ET AL. ***
      g_NDS = sqrt(20.448 * 4.*pi);    //! as N-D*-Sigma*_c
      gg = (2./3.) * pow(g_NDS,2.) / (16.0 * pi*pi * mRo*mRo); //  !WE INCLUDE AN OVERALL FACTOR OF 2/3 TO ACCOUNT FOR ISOSPIN!!!!
    }
    else if(dis==1){
      //!!**** A BLANK SPACE FOR OTHER SEPARATE MODES
      mRo = 0.7754;     //  !THE MASS OF THE RHO MESON
      mD = 1.232;       // !Mass OF THE DELTA BARYON
      //!!*** TERMINATE THE OUTPUT ***
      g_NDS = 0.0;      //! ZERO OUTPUT
      gg = pow(g_NDS,2.) / (16.0 * pi*pi * mRo*mRo);
      // !WE INCLUDE AN OVERALL FACTOR OF 1 TO ACCOUNT FOR ISOSPIN!!!!
    }

  if(y<=0. || y>=0.999){
    f_RhoDel = 0.0;
    return f_RhoDel;
  }     
  kT2 = kT*kT;
  SRoD = (kT2 + mRo*mRo)/y + (kT2 + mD*mD)/(1.0-y);
        
  if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SRoD));           //! monopole
  else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + SRoD)),2.);    //! dipole
  else if(typ==2)  FF = exp( (mN*mN - SRoD)/(L*L) );              //! expon
  else if(typ==3){
    t = ( -kT2 - (1.0-y)*(mRo*mRo - y*mN*mN) ) / y;
    FF = pow(((L*L - mD*mD) / (L*L - t)),2.);                      //! t-channel dip
  }
  else if(typ==4){
    sM = (kT2 + (2.0-y)*mRo*mRo)/(1.0-y) + (kT2 + (1.0-y)*mD*mD)/y + mN*mN;
    FF = (pow(L,4.) + pow(mD,4.))/(pow(L,4.) + pow(sM,2.));            //!DIPOLE -- s-channel Lambda exchange
  }


  //...TOPT with P[alpha] - p[alpha] derivative coupling
  P_k = (mRo*mRo + y*y*mN*mN + kT2)/2.0/y;
  P_p = (mD*mD + pow((1.0-y),2.)*mN*mN + kT2)/2.0/(1.0-y);
  pl_k = (mD*mD+kT2)*y/2.0/(1.0-y) + (mRo*mRo+kT2)*(1.0-y)/2.0/y + kT2;  

         sr = -4.0*mN*mD/3.0*(2.0*mD*mD+mN*mD+2.0*mN*mN)
	   -4.0*mN*mD/(3.0*mRo*mRo)*pow((P_k-pl_k),2.)
        -4.0/(3.0*mRo*mRo)*(mD*mD*P_k*P_k+mN*mN*pl_k*pl_k)
        +4.0*P_p/3.0*(2.0*mD*mD+4.0*mN*mD+mN*mN)
	   +4.0*P_p/(3.0*mRo*mRo)*pl_k*pl_k*(1.0-mN*mN/(mD*mD))
        -4.0*P_p*P_p*(1.0-2.0*P_k*pl_k/(3.0*mRo*mRo*mD*mD)
		       -P_p/(3.0*mD*mD));
          
	 if(kT>0.0 && (sr<0.0)){
	   cout << "CS2 -- ##### kT,y,sr =" << kT << ", "<< y << ", "<< sr << endl;
	   // should be stopped !!
	 }

	 ss = sr /pow(((1.0-y)*(SRoD - mN*mN)),2.) * FF*FF;

	 f_RhoDel = gg * (1.0-y) / y * ss;

	 return f_RhoDel;
}











// ******************************************************************************
/*
 * Function giving numerical value of f(y) for the SU(4) analogue of the pi-Delta
 * interaction, as usual y is IMF momentum fraction of the charmed BARYON.
*/
// ******************************************************************************

double fypiD(double y, double kT, double L, int typ, int dis){
  double ss,kT2,SpiD;
  double yR,fypiD,t,sM;
  double pi,mN,mD,mpi=0,g_DLcN,gg=0,FF=0;

  
  pi = 4*atan(1.);
  mN  = 0.93891897; //masses in GeV!!

  if (dis==0){
    //**** THE DISSOCIATION N ---> (pi^0 + DELTA^0) + (pi^- + DELTA^-)
    mD = 1.232;        //  !THE MASS OF THE DELTA BARYON
    mpi = 0.13957018;  // !Mass OF THE PION 
    //*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
    g_DLcN = sqrt (0.2237 * 4*pi);    //! as N-pi-Del from Holzenkamp et al.
    gg =  (2./3.) * pow(g_DLcN,2.) / (16. * pi*pi);   //!2/3 ISOSPIN FACTOR
  }
  else if (dis==1){
    //!!**** THE NULL DISSOCIATION
    mD = 1.232;         //  !THE MASS OF THE DELTA
    mpi = 0.13957018;   // !Mass OF THE PION
    //!!*** WE USE ZERO COUPLINGS TO RETURN A NULL OUTPUT
    g_DLcN = 0.;    
    gg = pow(g_DLcN,2.) / (16.0 * pi*pi);
  }

  if(y<=0.0 || y>=0.999){
    fypiD = 0.;
    return fypiD;
  }
  kT2 = kT*kT;
  SpiD = (kT2 + mpi*mpi)/y + (kT2 + mD*mD)/(1.-y);

  
  if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SpiD));         //! monopole
  else if (typ==1) FF = pow(((L*L + mN*mN) / (L*L + SpiD)),2.);  //! dipole
  else if (typ==2) FF = exp( (mN*mN - SpiD)/(L*L) );            //! expon
  else if(typ==3){ 
    t = (- kT2 - pow(mN,2.)*y*y) /(1.-y);
    FF = pow(((L*L - mpi*mpi) / (L*L - t)),2.);                  //! cov dip
  }
  else if(typ==4){                                              //! DIPOLE -- s-channel Lambda exchange
    sM = (kT2 + (1.+y)*mpi*mpi)/y + (kT2 + y*mD*mD)/(1.-y) + mN*mN;
    FF = (pow(L,4.) + pow(mD,4.))/(pow(L,4.) + pow(sM,2.));
  }

  yR = 1. - y;
  ss = ( kT2 + pow((mD-yR*mN),2.)) * pow(( kT2 + pow((mD+yR*mN),2.) ),2) / ( 6.*mD*mD*pow(yR,3.)) / pow(((1.-yR)*(SpiD - mN*mN)),2.) * FF*FF;

  fypiD =  gg/(mpi*mpi) * (1.0-yR) / yR * ss;

  return fypiD;

}





// ***************************************************************************
//        THE PION STRUCTURE FUNCTION IS TAKEN FROM THIS SUBROUTINE
/*
 *  Subroutine giving x-dependence of parametrization of LO pion valence
 *  and sea distribution functions xVpi, xSpi, valid for 0.25 < Q2 < 10^8 GeV^2,
 *  and 10^-5 < x < 1.
 *
 *  Gluck, Reya, Vogt: Z.Phys. C53 (1992) 651 (appendix 1).
 */
// ***************************************************************************
double GRV_xVpi(double x,double Q2){

  double  a,AA,D,Nv,alpha,as,AAs,Bs,Ds,E,Epr,beta;
  double  Q02,L,s;

  double xVpi = 0.0;
  double xSpi = 0.0;
 
  if(x<1.e-5) x=1.01e-5;
  if(Q2<0.25) Q2=0.2501;

  Q02 = 0.25;
  L = 0.232;
  if(Q2<Q02) return xVpi=0.;
  s = log( log(Q2/L/L) / log(Q02/L/L) );
  //...Valence distribution
  Nv = 0.519 + 0.180*s - 0.011*s*s;
  a = 0.499 - 0.027*s;
  AA = 0.381 - 0.419*s;
  D = 0.367 + 0.563*s;
  xVpi = 0.;
  
  if(x<1.){
    xVpi = Nv * pow(x,a) * (1.+AA*sqrt(x)) * pow((1.-x),D);
       }
  //C...Sea distribution (SU(3) symmetric)
  alpha = 0.55;
  as = 2.538 - 0.763*s;
  AAs = -0.748;
  Bs = 0.313 + 0.935*s;
  Ds = 3.359;
  E = 4.433 + 1.301*s;
  Epr = 9.30 - 0.887*s;
  beta = 0.56;
  xSpi = 0.;
  
  if(x<1.0){
    xSpi = pow(s,alpha) / pow((log(1./x)),as)
      * (1. + AAs*sqrt(x) + Bs*x) * pow((1.- x),Ds)
      * exp(-E + sqrt(Epr*pow(s,beta)*log(1./x)));
  }
  
  return xVpi;
}



double GRV_xSpi(double x,double Q2){

  double  a,AA,D,Nv,alpha,as,AAs,Bs,Ds,E,Epr,beta;
  double  Q02,L,s;

  double xVpi = 0.0;
  double xSpi = 0.0;
 
  if(x<1.e-5) x=1.01e-5;
  if(Q2<0.25) Q2=0.2501;

  Q02 = 0.25;
  L = 0.232;
  if(Q2<Q02) return xSpi=0.;
  s = log( log(Q2/L/L) / log(Q02/L/L) );

  //...Valence distribution
  Nv = 0.519 + 0.180*s - 0.011*s*s;
  a = 0.499 - 0.027*s;
  AA = 0.381 - 0.419*s;
  D = 0.367 + 0.563*s;
  xVpi = 0.;
  
  if(x<1.){
    xVpi = Nv * pow(x,a) * (1.+AA*sqrt(x)) * pow((1.-x),D);
       }
  //C...Sea distribution (SU(3) symmetric)
  alpha = 0.55;
  as = 2.538 - 0.763*s;
  AAs = -0.748;
  Bs = 0.313 + 0.935*s;
  Ds = 3.359;
  E = 4.433 + 1.301*s;
  Epr = 9.30 - 0.887*s;
  beta = 0.56;
  xSpi = 0.;
  
  if(x<1.0){
    xSpi = pow(s,alpha) / pow((log(1./x)),as)
      * (1. + AAs*sqrt(x) + Bs*x) * pow((1.- x),Ds)
      * exp(-E + sqrt(Epr*pow(s,beta)*log(1./x)));
  }
  
  return xSpi;
}








