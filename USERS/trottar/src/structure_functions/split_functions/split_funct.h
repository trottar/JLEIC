// ******************************************************************************
/* FUNCTION GIVING f(y) FOR N-PION VERTEX, WHERE
 *    y IS THE IMF MOMENTUM OF THE INTERMEDIATE MESON.
 *
 *  WRITTEN: W. Melnitchouk (1999)
 *  MODIFIED: T. HOBBS (2013)
 */
// ******************************************************************************
double fypiN(double y,double kT,double L,int typ,int dis){

  double ss,kT2,SpiN;
  double fypiN,t;
  double pi,mN,mpi=0,mP=0,g_piNN,gg=0,FF=0,sM;

  pi = 4*atan(1.0);
  mN  = 0.93891897;  //!masses in GeV!!
  //!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
  g_piNN = sqrt(14.40 * 4*pi);                 //! g_{pi NN}
  gg = 2.0 * pow(g_piNN,2.) / (16.0 * pi*pi);   //!2 ISOSPIN FACTOR

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
/*  Function giving numerical value of f(y) for N-Lambda-K vertex
 *  y is l.c. momentum fraction on baryon.
 *  Written: W. Melnitchouk (1999)
 *  Modified: T. Hobbs (2012)
 ****************************************************************************
 */
double fykL(double y,double kT,double L,int typ,int dis){

  double ss,ss0,ikT,kT2,kTmax,kTint,SkL;
  double fykL,t;
  double pi,mN,mk=0,mL=0,g_kLN,gg=0,FF=0,sM;

  pi = 4*atan(1.0);
  mN  = 0.93891897;  //!masses in GeV!!

  mL = 1.1157; // Mass of LAMBDA^0
  mk = 0.4937; // Mass of K^+
  //!!*** WE USE THE COUPLINGS INFERRED FROM MUELLER-GROELING ET AL. ***
  g_kLN = sqrt(15.56 * 4*pi);    //! Mueller-Groeling - PS
  //g_kLN = sqrt(14.40 * 4*pi)    //! g_{pi NN}
  //g_kLN = sqrt(9. * 4*pi)       //! Navarra I
  //g_kLN = sqrt(3.61 * 4*pi)     //! Navarra II
  gg = g_kLN*g_kLN / (16. * pi*pi);
    
  ss0 = 0.;
  ikT=1;
  kTmax = 10.;
  kTint = kTmax/1000.;

  for(kT=kTint;kT<=kTmax;kT++){
    kT2 = kT*kT;
    SkL = (kT2 + mk*mk)/(1.-y) + (kT2 + mL*mL)/y;

    if(typ==0)       FF = ((L*L + mN*mN) / (L*L + SkL));             //! monopole
    else if(typ==1)  FF = pow(((L*L + mN*mN) / (L*L + kT)),2.) ;     //! dipole
    else if(typ==2)  FF = exp( (mN*mN - kT)/(L*L) );                //! expon
    else if(typ==3) {
      t = (- kT2 - mN*mN*y*y) /(1.0-y);
      FF = pow(((L*L - mk*mk) / (L*L - t)),2.);                      //! cov dip
    }
    else if(typ==4){                       
      sM = (kT2 + (1.0+y)*mk*mk)/y +
	(kT2 + y*mL*mL)/(1.0-y) + mN*mN;
      FF = (pow(L,4.) + pow(mL,4.))/(pow(L,4.) + sM*sM); // ! DIPOLE -- s-channel Lambda exchange
    }
    ss = ( kT2 + pow((mL - y*mN),2.)) / y
      / pow(( (1.0-y)*(SkL - mN*mN) ),2.) * FF*FF * (2*kT);

    if(ikT/2*2!=ikT) ss0 = ss0 + 4*ss;
    else if(ikT/2*2==ikT) ss0 = ss0 + 2*ss;
    else ikT = ikT +1; 
  }

  fykL =  gg * (1.0-y) / y * (kTint/3) * ss ;
  return fykL;

} 
