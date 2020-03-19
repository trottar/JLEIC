/* for the JLEIC configuration
 *  TDISMC_collMC
 *  Created by K.Park and R.Trotta <trotta@cua.edu>
 */

#include <stdio.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TClass.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <TRandom.h>
#include <assert.h>
#include <TBuffer.h>
#include "cteq/cteqpdf.h"

using std::cout;
using std::endl;
using std::scientific;
using std::fixed;
using std::ios;


int NEvts = 0; // Defined in kinematics.input

// spectator proton either proton beam or deuteron beam
const double pSMax=  0.3;


//  Physical Constants
const double MProton   = 0.93827203;
const double MNeutron  = 0.93956536;
const double mElectron = 0.5110e-3;
const double mPim     = 0.13957018;
const double mPip     = 0.13957018;
const double mPi0     = 0.1349766;
const double mKp      = 0.493667;
const double MLambda  = 1.1156;
const double MDeut     = MProton+MNeutron-0.0022;
const double MBindA    = -0.008; // Average Binding Energy
const double MPSq      = MProton*MProton;
const double alphaQED  = 1./137.03;
const double pi        = acos(-1.0);

// Initialize Beam
double PBeam = 0;  // ion Beam momentum/Z (GeV/c), Defined in kinematics.input
double kBeam =  0.;  // Electron Beam Momentum, Defined in kinematics.input
//  electron and ion beam polarization
const double eBeamPol = 1.;
const double DBeamPol = 1.;


const  double PI     = 3.1415;
const  double D2R    = PI/(180.0);


const double ZBeam =   1.;
const double ABeam =   1.;           // 1 for proton beam, 2 for deuteron beam
//const double CrossingTheta = 0.050; // JLEIC angle
const double CrossingTheta = +0.025; // eRHIC angle
//0.050; // 0.010 for eRHIC/ePHENIX
const double CrossingPhi   = 0.0;   // pi/2.0 for eRHIC/ePHENIX

const double eBetaStarX = 0.10; // electron, ion beta at IP (m)
const double eBetaStarY = 0.02;
const double iBetaStarX = 0.10;
const double iBetaStarY = 0.02;

// Normalized emittance values (m*radian) nominal setup : Transverse
const double eEpsNX     = 54.e-6; 
const double eEpsNY     = 11.e-6;
// geometrical emittance iEpsX = M_i/P_i * iEpsNX
const double iEpsNX     = 0.35e-6;
const double iEpsNY     = 0.07e-6;

const double eDkOverk   = 7.1e-4; // Fractional energy spread  Normalized emittance values : Longitudinal
const double iDPoverP   = 3.0e-4;




// Global event-by-event Invariant Variable Structure
const double nanobarn = 1.0e-33; // cm2



// beam smearing function call
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt){
  double gamma = sqrt(pInc*pInc+mInc*mInc);
  double sig   = sqrt(NormEmit/(gamma*betaSt));
  return sig;
}




// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 6 ! t-dependent expon FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pipaulvillas(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.37229E-03;
    p1 = 0.69442E-02;
    p2 = -0.16933;
    p3 = 0.31326;
    p4 = 6.6281;
    p5 = -26.467;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.37229E-03;
    p1 = 0.69442E-02;
    p2 = -0.16933;
    p3 = 0.31326;
    p4 = 6.6281;
    p5 = -26.467;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.42299E-03;
    p1 = 0.33955E-01;
    p2 = -0.64850;
    p3 = 4.3910;
    p4 = -13.446;
    p5 = 15.984;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.89681E-03;
    p1 = 0.60363E-01;
    p2 = -0.88363;
    p3 = 4.7711;
    p4 = -11.958;
    p5 = 11.768;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ // Hall-A fixed target exp. limit */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 5 ! t-dependent expon FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pitexp(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.13151E-04;
    p1 = 0.24655E-02;
    p2 = -0.11011;
    p3 = 1.6664;
    p4 = -10.813;
    p5 = 25.662;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.21109E-03;
    p1 = 0.24422E-01;
    p2 = -0.65408;
    p3 = 6.1526;
    p4 = -25.771;
    p5 = 41.141;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.53755E-03;
    p1 = 0.42349E-01;
    p2 = -0.80465;
    p3 = 5.4029;
    p4 = -16.362;
    p5 = 19.186;
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.13141E-02;
    p1 = 0.89010E-01;
    p2 = -1.3003;
    p3 = 7.0225;
    p4 = -17.631;
    p5 = 17.393;
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 3 ! COV DIP FORM FACTOR
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2picov(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.30433e-04;
    p1 = 0.60118e-03;
    p2 = -0.30863e-01;
    p3 = 0.74473e-01;
    p4 = 4.2508;
    p5 = -28.398;
    //    if (x < 0.0555 || x > 0.083) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.43685e-03;
    p1 = 0.83395e-02;
    p2 = -0.20040;
    p3 = 0.36447;
    p4 = 7.8242;
    p5 = -31.101;
    //    if (x < 0.0555 || x > 0.16) xflag = 1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.60436e-03;
    p1 = 0.45603e-01;
    p2 = -0.86055;
    p3 = 5.6867;
    p4 = -16.812;
    p5 = 19.095;
    //     if (x < 0.0555 || x > 0.226) xflag =1;
    // Comment out for EIC
    /* if ( x > 0.3) xflag = 1; */
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.16068e-02;
    p1 = 0.10964;
    p2 = -1.5980;
    p3 = 8.6327;
    p4 = -21.714;
    p5 = 21.475;
    //     if (x < 0.0555 || x > 0.281) xflag = 1;
    if ( x > 0.3) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }

  // Comment out for EIC
  /* if( xflag == 1 || x < 0.0555 || x > 0.3){ */
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta (Timothy J. Hobbs)
// what is "th" = theta angle
//
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 2 !EXPONENTIAL FORM FACTOR (s-depdendent exp)
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double f2pi(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 3.656e-5;
    p1 = -0.000402;
    p2 = -0.008886;
    p3 = 0.07359;
    p4 = 1.079;
    p5 = -8.953;
    //if (x < 0.0555 || x > 0.083) xflag = 1;
    //    fk = -0.954 + 66.5*p -1632.4*p*p + 14573*p*p*p; // original value from Dasuni
    fk = -0.2473 + 14.79*p -283.55*p*p + 1766.9*p*p*p;   // fitting 3Var_k.f  x1000
  }

  if (p > 0.1 && p <= 0.2){
    p0 = 0.000287;
    p1 = 0.009397;
    p2 = -0.2632;
    p3 = 2.029;
    p4 = -5.878;
    p5 = 4.664;
    //if (x < 0.0555 || x > 0.16) xflag = 1;
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
    p0 = 0.0003662;
    p1 = 0.02811;
    p2 = -0.4566;
    p3 = 2.428;
    p4 = -5.107;
    p5 = 3.222;
    // if (x < 0.0555 || x > 0.226) xflag =1;
    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.0009412;
    p1 = 0.01366;
    p2 = -0.1744;
    p3 = 0.3864;
    p4 = 0.6615;
    p5 = -2.113;
    // if (x < 0.0555 || x > 0.281) xflag = 1;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }
  
  //  if( xflag == 1 || x < 0.0555 || x > 0.3){
  if( xflag == 1 ){      // EIC collider no limit
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
}

// subroutine to calculate the f2pi as function of xbj (Timothy J. Hobbs)
// using the ZEUS parameterization with f2p
double f2piZEUS(double x){
  double f2 = 0.0;
  double p0 = 0.1736;
  double p1 = 4.537;
  double p2 = -48.66;
  double p3 = 236.8;
  double p4 = -665.8;
  double p5 = 1094;
  double p6 = -973.9;
  double p7 = 363.2;

  // EIC collider no limit
  /* if( x < 0.0 || x > 0.6 )  */
    /* return 0.0;     */
  /* else{ */
  f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5) + p6*pow(x,6) + p7*pow(x,7);
  double f2temp = 0.361*f2;
  return f2temp;
  /* } */
}


// Testing K. Park  Nov. 08 2016
// subroutine to calculate the f2kp as function of recoiled nucleon momentum, xbj, theta
// This is user parametrization by fit the Wally's codes 3KVar_x.f() with integration of finite momentum range
// typ = 2 ! s-exp For factor
// dis = 0 ! Charge EXCHANGE
// FLAG = 0  --- THE Kaon CONTRIBUTION      | J = 0 + 1/2

double f2kp(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
      p0 = 0.83443E-07;
      p1 = -0.13408E-05;
      p2 = 0.22470E-03;
      p3 = -0.75778E-02;
      p4 = 0.85898E-01;
      p5 = -0.32083;
      //if (x < 0.0555 || x > 0.083) xflag = 1;
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p;
    //    cout << "--> momentum region 1: "  << endl;  
  }

  if (p > 0.1 && p <= 0.2){
      p0 = 0.32072E-05;
      p1 = 0.50358E-03;
      p2 = -0.76921E-02;
      p3 = 0.23637E-01;
      p4 = 0.87552E-01;
      p5 = -0.39142;
      //if (x < 0.0555 || x > 0.16) xflag = 1;
    fk = 0.464 -15.4*p + 126.5*p*p;
    //    cout << "--> momentum region 2: "  << endl;
  }
  
  if (p > 0.2 && p <= 0.3){
      p0 = 0.26275E-04;
      p1 = 0.30344E-02;
      p2 = -0.38997E-01;
      p3 = 0.15174;
      p4 = -0.16257;
      p5 = -0.80353E-01;
      //if (x < 0.0555 || x > 0.226) xflag =1;
     fk = -1.133 + 8.5354*p;
     //    cout << "--> momentum region 3: "  << endl;
  }

  if (p > 0.3 && p <= 0.5){
    p0 = 0.18012E-03;
    p1 = 0.18799E-01;
    p2 = -0.20900;
    p3 = 0.88713;
    p4 = -1.8633;
    p5 = 1.7046;
    // if (x < 0.0555 || x > 0.281) xflag = 1;
     fk = -1.345 + 9.47*p -7.91*p*p;
     //    cout << "--> momentum region 4: "  << endl;
  }


  if (p > 0.5 && p <= 1.){
    p0 = 0.12115E-02;
    p1 = 0.72634E-01;
    p2 = -0.60468;
    p3 = 2.4920;
    p4 = -6.1261;
    p5 = 6.3972;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    //    cout << "--> momentum region 5: "  << endl;
  }

  
  if (p > 1.0 && p <= 5.){
    p0 = 0.50862E-03;
    p1 = 0.30568E-01;
    p2 = -0.23605;
    p3 = 0.98762;
    p4 = -2.4529;
    p5 =  2.5330;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    //    cout << "--> momentum region 6: "  << endl;
  }



  if (p > 5.0 && p <= 100.){
    p0 = 0.20262E-06;
    p1 = 0.25285E-04;
    p2 = -0.45309E-03;
    p3 = 0.28817E-02;
    p4 = -0.80264E-02;
    p5 =  0.83156E-02;
    fk = -1.345 + 9.47*p -7.91*p*p; // same as p=0.3-0.5 GeV/c 
    //    cout << "--> momentum region 7: "  << endl;
  }

  

  
  if (p < 0.05){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
    //    cout << "--> momentum region 0: "  << endl;
  }
  

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }
  
  //  if( xflag == 1 || x < 0.0555 || x > 0.3){
  if( xflag == 1){
    return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
 }


// typ = 3 ! t-exp For factor
// dis = 0 ! Charge EXCHANGE
// FLAG = 0  --- THE Kaon CONTRIBUTION      | J = 0 + 1/2

double f2kptmono(double p, double x, double th){
  
  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    p0 = 0.20488E-05;
    p1 = 0.63463E-03;
    p2 = -0.27479E-01; 
    p3 = 0.41698;
    p4 = -2.7294;
    p5 = 6.5210;
    //if (x < 0.0555 || x > 0.083) xflag = 1;
    if (x > 0.08) xflag = 1;
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573.*p*p*p; 
  }

  if (p > 0.1 && p <= 0.2){
      p0 = 0.41053E-04;
      p1 = 0.55437E-02;
      p2 = -0.14278;
      p3 = 1.3096;
      p4 = -5.3818;
      p5 = 8.4883;
      //if (x < 0.0555 || x > 0.16) xflag = 1;
      if ( x > 0.2) xflag = 1;
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  
  if (p > 0.2 && p <= 0.3){
      p0 = 0.16415E-03;
      p1 = 0.13833E-01;
      p2 = -0.25650;
      p3 = 1.6955;
      p4 = -5.0825;
      p5 = 5.9353;
      //if (x < 0.0555 || x > 0.226) xflag =1;
      if ( x > 0.226) xflag =1;
     fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
      p0 = 0.69820E-03;
      p1 = 0.48810E-01;
      p2 = -0.70672;
      p3 = 3.8132;
      p4 = -9.6119;
      p5 = 9.5429;
      //if (x < 0.0555 || x > 0.281) xflag = 1;
      if ( x > 0.28) xflag = 1;
     fk = -1.345 + 9.47*p -7.91*p*p;
  }

  if (p < 0.05 || p > 0.5){
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;
  }

  if (th < 1.8 || th > 74){
    fth = 0.0;}
  else{
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; }
  
  if( xflag == 1 || x < 0.001 || x > 0.3){
      //if( xflag == 1 ){ // testing 
      return 0.0;}
  else{
    double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);  
    double f2temp = f2*fk*fth;
    return f2temp;
  }
  
  
 }


// subroutine to calculate the f2p as a function xbj
double f2p( double x ){
  double f2 = 0.0;
  double p0 = 0.1736;
  double p1 = 4.537;
  double p2 = -48.66;
  double p3 = 236.8;
  double p4 = -665.8;
  double p5 = 1094;
  double p6 = -973.9;
  double p7 = 363.2;

  f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5) + p6*pow(x,6) + p7*pow(x,7);
  
  return f2;
}

// pimake with smearing of fermi motion
TVector3 TVector3::Unit() const {

  ifstream in;
  in.open(Form("moment_ld2b.dat"));

  double degr=0.01745329252;
  double fgev=0.1973;

  const int nfermi = 201;
  double ppfermi[nfermi] = {0.0};
  double pfdis[nfermi]= {0.0};
  double pfdis1[nfermi] = {0.0};
  double pfdis2[nfermi]= {0.0};
  int i = 0;
  
  while(1){
    in>>ppfermi[i]>>pfdis[i]>>pfdis1[i]>>pfdis2[i];
    //cout << ppfermi[i] << "\t" << pfdis[i] << endl;
    ppfermi[i] = ppfermi[i]*fgev;
    pfdis[i] = pfdis[i]*ppfermi[i]*ppfermi[i];
    i++;
    if(i==nfermi) break;
  }
  in.close();
 

  for (int j = 1; j < nfermi; j++){
    pfdis[j] = pfdis[j] +  pfdis[j-1];
    //cout << ppfermi[j] << "\t" << pfdis[j] << endl;
  }

  if (pfdis[nfermi-1] > 0.0){
    for (int j = 0; j < nfermi; j++){
      pfdis[j] = pfdis[j]/pfdis[nfermi-1];
    }
  }
  else
    cout << "error in reading fermi file" << endl;
  
  double xran = gRandom->Uniform(1.);

  double fermi = -1.0;
  for (int i = 0; i < nfermi; i++){
    if (pfdis[i] == xran){
      fermi = ppfermi[i];
    }
    else if (xran > pfdis[i] && xran < pfdis[i+1]){
      double numerator   = ( ppfermi[i+1] -  ppfermi[i]);
      double denominator = (  pfdis[i+1] -   pfdis[i]);
      double slope = numerator/denominator;
      fermi = ppfermi[i] + (xran - pfdis[i])*slope;
     
    }
    
  }

  // Random number generation for xran1 and xran2 between [0.0 : 1.0]
  double xran1 = gRandom->Uniform(1.);
  double fth = acos(2*xran1 - 1.0);

  double xran2 = gRandom->Uniform(1.);
  double fphi = 2*180.0*xran2*degr;
  
  
  double fermix = fermi*sin(fth)*cos(fphi);
  double fermiy = fermi*sin(fth)*sin(fphi);
  double fermiz = fermi*cos(fth);
  

  // For debugging purpose: OK
  //  cout << "inside of Unit function:: xran1= " << xran1 << ", xran2= "
  //        << xran2 <<  " fermi= " << fermi <<  ", fth= " << fth <<  ", fphi= " << fphi << endl;
 
  TVector3 PiMake(fermix, fermiy, fermiz);

  return PiMake;
 
}

// Use CTEQ6 parameterization
cteq_pdf_t *__dis_pdf;


// ERROR: This or it being called below in F2N, the cteq_pdf_alloc_id
// function is opening and closing a file for every event, using up
// resources.
void initcteqpdf(){
  __dis_pdf = cteq_pdf_alloc_id(400); // mode 400 = cteq6.6?

  assert(__dis_pdf);
}

double F2N(double x, double Q2, int nucl){
  // nucl =2; // for neutron;

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  /* initcteqpdf(); */
  
  double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

  double quv = qu-qubar;
  double qdv = qd-qdbar;

  double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;

  if( nucl == 1 ){
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
  }
  if( nucl == 2){
    F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
  }
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
  return F2;

}

// DIS cross-section from CTEQ parameters
//double dissigma( double ebeam, double th, double eprime, int nucl ){ // ** fixed target exp.
/*
  double Q2 = 2.0*eprime*ebeam*(1.0-cos(th));
  double nu = ebeam-eprime;
  double Mp = 0.938;

  double x = Q2/(2.0*Mp*nu);
  double y = nu/ebeam;
*/

// for the JLEIC collider : it needs an input as x, y_D, q2, nu, nucl instead of ebeam, theta_e, eprime 
double cdissigma( double k_x, double k_y, double k_q2, double k_nu, double k_ep, int nucl ){ // ** collider exp.
  // Return in nb
  double Q2 = k_q2;
  double nu = k_nu;
  double Mp = 0.938;
  double eprime = k_ep;
  double x = k_x;
  double y = k_y;

  // for debugging purpose
  /*
    cout << "x = " << x << ", y= " << y << ", nu= " << nu << ", Q2= " << Q2 << endl;
  */

  if( ! (0.0 < x && x < 1.0 && 0.0 < y && y < 1.0) ){
    //printf("WARNING %s line %d  x = %f, y = %f -> eprime = %f GeV   th = %f deg  ebeam = %f GeV\n", __FILE__,
    //__LINE__, x, y, eprime, th*180/3.14159, ebeam );
    //exit(1);
    return 0.0;;

  }

  double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

  double quv = qu-qubar;
  double qdv = qd-qdbar;

  double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;

  if( nucl == 1 ){
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
  }
  if( nucl == 2){
    F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
  }
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
  double F1 = F2/(2.0*x);

  // From PDG
  double ds_dxdy= 4.0*3.14159*((1.0-y-pow(x*y*Mp,2.0)/Q2)*F2+y*y*x*F1)/(x*y*Q2*137.0*137.0);

  // In GeV^-2
  double ds_dOmega_dE = ds_dxdy*eprime/(2.0*3.14159*Mp*nu);

  // From Zeus paper
  double ds_dxdQ2 = (2*(3.14159)/(x*Q2*Q2*(137)*(137)))*((1+((1-y)*(1-y)))*(F2)-(1+((1-y)*(1-y)))*(x*F1));
    
  /* return ds_dOmega_dE*0.197*0.197*1e7; // GeV2 -> nb */

  return ds_dxdQ2*0.197*0.197*1e7; // GeV2 -> nb 
}

// Call specific cross-section function w.r.t target nucleon
/*
  double dissigma_p(double eb, double th, double ep){
  return dissigma( eb, th, ep, 1);
  }
  double dissigma_n(double eb, double th, double ep){
  return dissigma( eb, th, ep, 2 );
  }
*/

// *** for the JLEIC collider : it needs an input as x, y_D, q2, nu instead of ebeam, theta_e, eprime 
double cdissigma_p(double x, double y_D, double q2, double nu, double eprime){
  return cdissigma( x, y_D, q2, nu, eprime, 1);
}
double cdissigma_n(double x, double y_D, double q2, double nu, double eprime){
  return cdissigma( x, y_D/2, q2, nu, eprime, 2);
}

/*
 *   Call CTEQ PDF function listed 
 */

/*    some technical constatnts  */
static const double xpow = 0.3;
static const double onep = 1.00001;
static const unsigned int maxval = 4;

static
double __cteq_pdf_as(unsigned int ord, unsigned int nf, double q, double lmd)
{
  double b0 = 11.0 - 2.0/3.0*nf;
  double t = log(q/lmd);
  double as = 1.0/(b0*t);
  
  /*  it is a leading order evolution  */ 
  if(ord <= 1) return as;
  
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double b1 = 51.0 - 19.0/3.0*nf;
  return as*(1.0-b1/b0*as*log(2.0*t));
}

static 
double __cteq_pdf_astolmd(unsigned int ord, unsigned int nf, double as, double q)
{
  double b0 = 11.0 - 2.0/3.0*nf;
  double t = 1.0/(b0*as);
  
  /*  it is a leading order evolution  */ 
  if(ord <= 1) return q*exp(-t);
  
  /*  at NLO or higer order level returns with the NLO alpha_s  */
  double as0, as1, ot, lt, br = (51.0 - 19.0/3.0*nf)/(b0*b0);
  
  do {
    lt = log(2.0*t)/t;
    ot = t;
	
    as0 = (1.0 - br*lt)/(b0*t);
    as1 = (-1.0 - br*(1.0/t-2.0*lt))/(b0*t*t);
    t += (as - as0)/as1;
  } while(fabs(ot-t)/ot > 1e-5);
  
  return q*exp(-t);
}

static
void __cteq_pdf_setlmd(cteq_pdf_t *pdf)
{
  double as;
  unsigned int nf;
  
  for(nf = pdf->nf+1; nf <= 6; ++nf) {
    as = __cteq_pdf_as(pdf->order, nf-1, pdf->mass[nf], pdf->lambda[nf-1]);
    pdf->lambda[nf] = __cteq_pdf_astolmd(pdf->order, nf, as, pdf->mass[nf]);
  }
  
  /*   Under the charm mass every quark is considered as massless.  */
  for(nf = pdf->nf-1; nf > 2; --nf) {
    as = __cteq_pdf_as(pdf->order, nf+1, pdf->mass[nf+1], pdf->lambda[nf+1]);
    pdf->lambda[nf] = __cteq_pdf_astolmd(pdf->order, nf, as, pdf->mass[nf+1]);
  }
}

static
double __cteq_pdf_polint4f(double *xa, double *ya, double x)
{
  double c1, d1, d2, c2, d3, h1, h2, h3, h4, c3, cc1, cd1, 
    cd2, cc2, dd1, dc1, den;
  
  /* Function Body */
  h1 = xa[0] - x; h2 = xa[1] - x;
  h3 = xa[2] - x; h4 = xa[3] - x;
  
  den = (ya[1] - ya[0])/(h1 - h2);
  d1 = h2*den; c1 = h1*den;
  
  den = (ya[2] - ya[1])/(h2 - h3);
  d2 = h3*den; c2 = h2*den;
  
  den = (ya[3] - ya[2])/(h3 - h4);
  d3 = h4*den; c3 = h3*den;
  
  den = (c2 - d1)/(h1 - h3);
  cd1 = h3*den; cc1 = h1*den;
  
  den = (c3 - d2)/(h2 - h4);
  cd2 = h4*den; cc2 = h2*den;
  
  den = (cc2 - cd1)/(h1 - h4);
  dd1 = h4*den; dc1 = h1*den;
  
  if(h3 + h4 < 0.0) return ya[3] + d3 + cd2 + dd1;
  if(h2 + h3 < 0.0) return ya[2] + d2 + cd1 + dc1;
  if(h1 + h2 < 0.0) return ya[1] + c2 + cd1 + dc1;
  
  return ya[0] + c1 + cc1 + dc1;
} 

static 
double __cteq_pdf_pardis(const cteq_pdf_t *pdf, int iprtn, double x, double q)
{
  int jm, jx, jlx = -1, ju, jq, jlq = -1, j1, ip, jtmp, it;
  unsigned int nx = pdf->nx, nq = pdf->nt;
  double tt, ss, fvec[4];
  
  if(iprtn != 0)
    if(q <= pdf->mass[abs(iprtn)]) 
      return 0.0;
  
  tt = log(log(q/pdf->lambda[pdf->nf]));
  ss = pow(x, xpow);
  
  ju = nx + 1;
  while(ju - jlx > 1) {
    jm = (ju + jlx)/2;
    if(x >= pdf->xv[jm]) jlx = jm;
    else ju = jm;
  }
  
  if(jlx <= -1) {
    fprintf(stderr, "cteq-pdf: Severe error x <= 0! x = %g \n", x);
    exit(-1);
  } else if(jlx == 0) jx = 0;
  else if(jlx <= (int) nx-2) jx = jlx-1;
  else if(jlx == (int) nx-1 || x < onep) jx = jlx-2;
  else {
    fprintf(stderr, "cteq-pdf: Severe error: x > 1!  x = %g", x);
    exit(-1);
  }
  
  ju = nq + 1;
  while(ju - jlq > 1) {
    jm = (ju + jlq)/2;
    if (tt >= pdf->tv[jm]) jlq = jm;
    else ju = jm;
  }
  
  if(jlq <= 0) jq = 0;
  else if(jlq <= (int) nq-2) jq = jlq-1;
  else jq = nq-3;
  
  /* get the pdf function values at the lattice points... */
  ip = (iprtn > (int) pdf->mxval ? -iprtn : iprtn);
  jtmp = ((ip + pdf->nfmx)*(nq+1) + (jq-1))*(nx+1) + jx + 1;
  
  if(jx == 0) 
    {
      double fij[4];	
      for(it = 0; it < 4; ++it) {
	j1 = jtmp + (it+1)*(nx+1);
	  
	fij[0] = 0.0;
	fij[1] = (pdf->upd[j1  ])*(pdf->xv[1])*(pdf->xv[1]);
	fij[2] = (pdf->upd[j1+1])*(pdf->xv[2])*(pdf->xv[2]);
	fij[3] = (pdf->upd[j1+2])*(pdf->xv[3])*(pdf->xv[3]);
	  
	/*  Use Polint which allows x to be anywhere w.r.t. the grid */
	fvec[it] = __cteq_pdf_polint4f(pdf->xvpow, fij, ss)/(x*x);
      }	
    } 
  else if((unsigned)jlx == nx-1) 
    {
      for(it = 0; it < 4; ++it)
	fvec[it] = __cteq_pdf_polint4f(pdf->xvpow+nx-3, pdf->upd+jtmp+(it+1)*(nx+1)-1, ss);
    } 
  else 
    {
      double *svec = pdf->xvpow + jx-1;
      double s12 = svec[1]-svec[2], s13 = svec[1]-svec[3], s23 = svec[2]-svec[3],
	s24 = svec[2]-svec[4], s34 = svec[3]-svec[4], sy2 = ss-svec[2], sy3 = ss-svec[3];
	
      double const1 = s13/s23, const2 = s12/s23, const3 = s34/s23, const4 = s24/s23;
      double s1213 = s12 + s13, s2434 = s24 + s34;
      double sdet = s12*s34 - s1213*s2434;
      double tmp = sy2*sy3/sdet; 
      double const5 = (s34*sy2 - s2434*sy3)*tmp/s12, const6 = (s1213*sy2 - s12*sy3)*tmp/s34;
	
      for(it = 0; it < 4; ++it) {
	j1 = jtmp + (it+1)*(nx+1);
	  
	double sf2 = pdf->upd[j1], sf3 = pdf->upd[j1+1];
	double g1 = sf2*const1 - sf3*const2, g4 = sf3*const4 - sf2*const3;
	fvec[it] = (const5*(pdf->upd[j1-1]-g1) + const6*(pdf->upd[j1+2]-g4) + sf2*sy3-sf3*sy2)/s23;
      }
    }
  
  /*   interpolate in t... */
  if(jlq <= 0) return __cteq_pdf_polint4f(pdf->tv, fvec, tt);
  if((unsigned)jlq >= nq-1) return __cteq_pdf_polint4f(pdf->tv+nq-3, fvec, tt);
  
  double *tvec = pdf->tv + jq-1;
  double t12 = tvec[1]-tvec[2], t13 = tvec[1]-tvec[3], t23 = tvec[2]-tvec[3],
    t24 = tvec[2]-tvec[4], t34 = tvec[3]-tvec[4], ty2 = tt-tvec[2], ty3 = tt-tvec[3];
  
  double tmp1 = t12 + t13, tmp2 = t24 + t34;
  double tdet = t12*t34 - tmp1*tmp2;

  double g1 = (fvec[1]*t13 - fvec[2]*t12)/t23, g4 = (fvec[2]*t24 - fvec[1]*t34)/t23;
  double h00 = (t34*ty2 - tmp2*ty3)*(fvec[0]-g1)/t12 + (tmp1*ty2 - t12*ty3)*(fvec[3]-g4)/t34;
  
  return (h00*ty2*ty3/tdet + fvec[1]*ty3 - fvec[2]*ty2)/t23;
}
  
static 
cteq_pdf_t * __cteq_pdf_alloc_read_tbl(FILE *file)
{
  unsigned int i;
  
  /*   Allocate the pdf set   */
  //cteq_pdf_t *pdf =  static_cast<cteq_pdf_t*> (malloc(sizeof(cteq_pdf_t)));
  cteq_pdf_t *pdf =  (cteq_pdf_t*) (malloc(sizeof(cteq_pdf_t)));
  if(!pdf) return 0;
  
  /*   Reading the first two line    */
  int del;
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading ord, nfl  */
  float __ord, __nf; 
  //fscanf(file, "%f%f", &__ord, &__nf);
  if (fscanf(file, "%f%f", &__ord, &__nf) == 1) {
    fprintf(file, "%f%f", __ord, __nf);
  } else {
    fprintf(file,"WARNING...\nFailed to read %f%f\n", __ord, __nf);
  }
  pdf->order = (unsigned int) __ord;
  pdf->nf = (unsigned int) __nf;
  pdf->mxval = 2;

  /*   Reading ord, nfl, al mass[1-6]   */
  pdf->mass[0] = 0.0;     /*   gluon mass  */
  fscanf(file, "%lf%lf%lf%lf%lf%lf%lf", pdf->lambda + pdf->nf, pdf->mass+1,pdf->mass+2,pdf->mass+3,pdf->mass+4,pdf->mass+5,pdf->mass+6);

  /*   Calculating the lambda values at the thresholds  */
  __cteq_pdf_setlmd(pdf);
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading nx, nt, nfmx   */
  fscanf(file, "%u%u%u", &(pdf->nx), &(pdf->nt), &(pdf->nfmx));
  
  /*   Allocating memory for the grid   */
  unsigned int nupd = (pdf->nx+1)*(pdf->nt+1)*(pdf->nfmx+1+pdf->mxval);

  if(!(pdf->xv = (double *) malloc((pdf->nx + 1)*sizeof(double)))) goto label_free_pdf;
  if(!(pdf->xvpow = (double *) malloc((pdf->nx +1)*sizeof(double)))) goto label_free_xv;
  if(!(pdf->tv = (double *) malloc((pdf->nt + 1)*sizeof(double)))) goto label_free_xvpow;
  if(!(pdf->upd = (double *) malloc(nupd*sizeof(double)))) goto label_free_tv;
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading qini nad qmax   */
  fscanf(file, "%le%le", &(pdf->qini), &(pdf->qmax));
  
  /*   Reading q values and converting to t  */
  double q;
  for(i = 0; i <= pdf->nt; i++) {
    fscanf(file, "%le", &q);
    pdf->tv[i] = log(log(q/(pdf->lambda[pdf->nf])));
  }

  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading xmin   */
  fscanf(file, "%le", &(pdf->xmin));

  /*   Reading x values and calculating xvpow values  */
  double x;
  for(i = 0; i <= pdf->nx; i++) {
    fscanf(file, "%le", &x);
    pdf->xv[i] = x;
    pdf->xvpow[i] = pow(x, xpow);
  }
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading the grid   */
  for(i = 0; i < nupd; i++)
    fscanf(file, "%le", pdf->upd+i);
  
  return pdf;

  /*  garbage collection   */
 label_free_tv: free(pdf->tv);
 label_free_xvpow: free(pdf->xvpow);
 label_free_xv: free(pdf->xv);
 label_free_pdf: free(pdf);
  
  return 0;
}

static 
cteq_pdf_t * __cteq_pdf_alloc_read_pds(FILE *file)
{
  unsigned int i;
  
  /*   Allocate the pdf set   */
  cteq_pdf_t *pdf =  (cteq_pdf_t*) (malloc(sizeof(cteq_pdf_t)));

  if(!pdf) return 0;
  
  /*   Reading the first two line    */
  int del;
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading ord, nfl   */
  float __ord, __nf;
  fscanf(file, "%f%f", &__ord, &__nf);
  
  pdf->order = (unsigned int) __ord;
  pdf->nf = (unsigned int) __nf;
  
  /*   Reading ord, nfl, al mass[1-6]   */
  pdf->mass[0] = 0.0;     /*   gluon mass  */
  fscanf(file, "%lf%lf%lf%lf%lf%lf%lf", pdf->lambda + pdf->nf, 
	 pdf->mass+1,pdf->mass+2,pdf->mass+3,pdf->mass+4,pdf->mass+5,pdf->mass+6);
  
  /*   Calculating the lambda values at the thresholds  */
  __cteq_pdf_setlmd(pdf);
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading nfmx and mxval   */
  fscanf(file, "%d%d%d%u%u", &del, &del, &del, &(pdf->nfmx), &(pdf->mxval));
  if(pdf->mxval > maxval) pdf->mxval = 3;
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');

  /*   Reading nx, n and ng   */
  unsigned int ng;
  fscanf(file, "%u%u%d%u", &(pdf->nx), &(pdf->nt), &del, &ng);
 
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Skipping the next ng lines and the comment starts QINI, QMAX,...,   */
  for(i = 0; i <= ng; i++) 
    do del = fgetc(file); while(del != '\n');
  
  /*   Allocating memory for the grid   */
  unsigned int nupd = (pdf->nx+1)*(pdf->nt+1)*(pdf->nfmx+1+pdf->mxval);

  if(!(pdf->xv = (double *) malloc((pdf->nx + 1)*sizeof(double)))) goto label_free_pdf;
  if(!(pdf->xvpow = (double *) malloc((pdf->nx +1)*sizeof(double)))) goto label_free_xv;
  if(!(pdf->tv = (double *) malloc((pdf->nt + 1)*sizeof(double)))) goto label_free_xvpow;
  if(!(pdf->upd = (double *) malloc(nupd*sizeof(double)))) goto label_free_tv;
    
  /*   Reading qini nad qmax   */
  fscanf(file, "%le%le", &(pdf->qini), &(pdf->qmax));
  
  /*   Reading t values */
  double tmp;
  for(i = 0; i <= pdf->nt; i++) 
    fscanf(file, "%le%le", &tmp, &(pdf->tv[i])); 
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading xmin   */
  fscanf(file, "%le%le", &(pdf->xmin), &tmp);
  
  /*   Reading x values and calculating xvpow values  */
  pdf->xv[i] = pdf->xvpow[i] = 0.0;
  
  for(i = 1; i <= pdf->nx; i++) {
    fscanf(file, "%le", &(pdf->xv[i]));
    pdf->xvpow[i] = pow(pdf->xv[i], xpow);
  }
  
  /*   Positioning and reading the comment line  */ 
  do del = fgetc(file); while(del != '\n');
  do del = fgetc(file); while(del != '\n');
  
  /*   Reading the grid   */
  for(i = 0; i < nupd; i++)
    fscanf(file, "%le", pdf->upd+i);
  
  return pdf;
  
  /*  garbage collection   */
 label_free_tv: free(pdf->tv);
 label_free_xvpow: free(pdf->xvpow);
 label_free_xv: free(pdf->xv);
 label_free_pdf: free(pdf);
  
  return 0;
}

/*    Exported functions    */
cteq_pdf_t * cteq_pdf_alloc(const cteq_pdfset_t *pdfset)
{
  /*   return value  */
  cteq_pdf_t *pdf = 0;
  
  /*   we need the full filename with the absolute path */
  size_t len = strlen(pdfset->path) + strlen(pdfset->filename);
  //char *filename = static_cast<char*>(malloc((len+1)*sizeof(char)));
  char *filename = (char*)(malloc((len+1)*sizeof(char)));
  
  if(!filename) {
    fprintf(stderr, "cteq_pdf: Unable to determine the full filename of the table file!\n");
    fprintf(stderr, "    path     : %s\n", pdfset->path);
    fprintf(stderr, "    filename : %s\n", pdfset->filename);
    return pdf;
  }

  strcpy(filename, pdfset->path);
  strcat(filename, pdfset->filename);
  
  /*   Openning the table file  */
  FILE *file = fopen(filename, "r");
  
  if(!file) {
    // Try top open locally stored file as a backup

    file = fopen(pdfset->filename, "r");

    if( !file ){
      fprintf(stderr, "cteq_pdf: Unable to open file %s or local file %s\n", filename, pdfset->filename);

      free(filename);
      return pdf;
    }
  }
  
  /*    Creating the pdf set   */
  switch(pdfset->itbl) {
  case 1: pdf = __cteq_pdf_alloc_read_tbl(file); break;
  case 2: pdf = __cteq_pdf_alloc_read_pds(file); break;
  default: pdf = 0;
  }
  
  /*    Closing the table file  */
  fclose(file);
  free(filename);
  
  return pdf;
}

cteq_pdf_t * cteq_pdf_alloc_name(const char *name){
  /*   First we try to find the pdfset in the database. */
  const cteq_pdfset_t *pdfset = cteq_pdfset_find(cteq_pdfset_database, name);
  if(pdfset) return cteq_pdf_alloc(pdfset);
  
  /*   Otherwise we try to open it as a table file. */
  /*   return value  */
  cteq_pdf_t *pdf = 0;
  
  /*   The table files must have .tbl or .pds extension. 
       Unfortunately there is no infor about the type of the table in the file,
       so we have to use the extension.   */
  unsigned int itbl = 0;
  const char *ext = name;
  while(*ext != '\0') ++ext; 
  
  if(strcmp(ext-4, ".tbl") == 0) itbl = 1;
  else if(strcmp(ext-4, ".pds") == 0) itbl = 2;
  else {
    fprintf(stderr, "cteq_pdf: Unable to identify the type of the table. Please use extension .tbl or .pds!\n");
    return pdf;
  }
  
  /*   Openning the table file  */
  FILE *file = fopen(name, "r");
  
  if(!file) {
    fprintf(stderr, "cteq_pdf: Unable to open file %s\n", name);
    return pdf;
  }
  
  /*    Creating the pdf set   */
  switch(itbl) {
  case 1: pdf = __cteq_pdf_alloc_read_tbl(file); break;
  case 2: pdf = __cteq_pdf_alloc_read_pds(file); break;
  default: pdf = 0;
  }
  
  /*    Closing the table file  */
  fclose(file);
  
  return pdf;
}

cteq_pdf_t * cteq_pdf_alloc_id(int id){
  /*   First we try to find the pdfset in the database. */
  const cteq_pdfset_t *pdfset = cteq_pdfset_find_id(cteq_pdfset_database, id);
  if(pdfset) return cteq_pdf_alloc(pdfset);
  
  return 0;
}

void cteq_pdf_free(cteq_pdf_t *pdf) 
{
  if(pdf == 0) return;
  
  if(pdf->xv) free(pdf->xv);
  if(pdf->tv) free(pdf->tv);
  if(pdf->xvpow) free(pdf->xvpow);
  if(pdf->upd) free(pdf->upd);
  
  free(pdf);
}

double cteq_pdf_evolvepdf(const cteq_pdf_t *pdf, int iprtn, double x, double q)
{
  if(abs(iprtn) > pdf->nfmx) return 0.0;  
  double ff = __cteq_pdf_pardis(pdf, iprtn, x, q);
  //  return ff <= 0.0 ? 0.0 : ff;
  if(ff <=0.) ff =0.;
  return ff;
}


double cteq_pdf_evolveas(const cteq_pdf_t *pdf, double q)
{
  unsigned int nf = 6;
  while(q < pdf->mass[nf] && nf > 3) --nf;
  return __cteq_pdf_as(pdf->order, nf, q, pdf->lambda[nf]);
}

/**   Returns the evolution order of the pdf set  */
unsigned int cteq_pdf_orderpdf(const cteq_pdf_t *pdf) {
  return pdf->order;
}

/**   Returns the evolution order of the \f$\alpha_s\f$ */
unsigned int cteq_pdf_orderas(const cteq_pdf_t *pdf) {
  return pdf->order;
}

/**   Returns the number of the active flavours  */
unsigned int cteq_pdf_nfmax(const cteq_pdf_t *pdf) {
  return pdf->nfmx;
}

/**   Returns the fitting scale  */
double cteq_pdf_scale(const cteq_pdf_t *pdf) {
  return pdf->qini;
}

/**   Returns the alphas at the fitting scale  */
double cteq_pdf_alfas(const cteq_pdf_t *pdf) {
  return cteq_pdf_evolveas(pdf, pdf->qini);
}

/**   Returns the masses  */
double cteq_pdf_mass(const cteq_pdf_t *pdf, int i) {
  return pdf->mass[abs(i)];
}

/**   Returns the flavour threshold in the evolution */
double cteq_pdf_threshold(const cteq_pdf_t *pdf, unsigned int nf) {
  return (nf < 4 ? 0.0 : pdf->mass[nf]);
}

//Add by Jixie: in case the makefile does not define CTEQ_TBL_PATH
//use this definition
#ifndef CTEQ_TBL_PATH
#define  CTEQ_TBL_PATH "cteq-tbls"
#endif

#define CTEQ6STD_TBL_PATH    CTEQ_TBL_PATH"/cteq6std/"
#define CTEQ6_TBL_PATH       CTEQ_TBL_PATH"/cteq6/"
#define CTEQ65S_TBL_PATH     CTEQ_TBL_PATH"/ctq65s/"
#define CTEQ65C_TBL_PATH     CTEQ_TBL_PATH"/ctq65c/"
#define CTEQ6M_TBL_PATH      CTEQ_TBL_PATH"/cteq6m/"
#define CTEQ61_TBL_PATH      CTEQ_TBL_PATH"/cteq61/"
#define CTEQ65_PDS_TBL_PATH  CTEQ_TBL_PATH"/ctq65-pds/"
#define CTEQ66A_TBL_PATH     CTEQ_TBL_PATH"/ctq66a/"
#define CTEQ66C_TBL_PATH     CTEQ_TBL_PATH"/ctq66c/"
#define CTEQ66M_TBL_PATH     CTEQ_TBL_PATH"/ctq66m/"

static const cteq_pdfset_t __cteq_pdfset_database[] = { 
  {400, "CTEQ66.00",  "description", CTEQ66M_TBL_PATH, "ctq66.00.pds", 2},
  
  /*  End of the list */
  {0,0,0,0,0,2}
};

const cteq_pdfset_t *cteq_pdfset_database = __cteq_pdfset_database;

const cteq_pdfset_t * 
cteq_pdfset_find(const cteq_pdfset_t *pdflist, const char *name)
{
  while(pdflist->name) {
    if(strcmp(name, pdflist->name) == 0) return pdflist;
    ++pdflist;
  }

  return 0;
}

const cteq_pdfset_t * 
cteq_pdfset_find_id(const cteq_pdfset_t *pdflist, int id)
{
  while(pdflist->name) {
    if(id == pdflist->id) return pdflist;
    ++pdflist;
  }
  
  return 0;
}

//Weightless uncertainty (for equal dist. of points)
double uncerWL(double val,int numData){

  double mean, stanDev = 0.0, numWeight = 0.0;
	
  mean = val/numData;
  stanDev += pow(val - mean, 2);
  //printf("Uncertainty is +/- %e for F2K \n", sqrt(stanDev/numData));
  /*
    for(int i=0; i<sizeof((int)weight[0]);i++){
    double numWeight = weight[i];
    return numWeight;
    }//*/
  return sqrt(stanDev/numData);
}

//Weighted uncertainty
double uncerW(double val,int numData,double weight[]){

  double mean, stanDev = 0.0, numWeight = 0.0;
	
  mean = val/numData;
  stanDev += pow(val - mean, 2);
  //printf("Uncertainty is +/- %e for F2K \n", sqrt(stanDev/numData));
  for(int i=0; i<(int)sizeof(weight[0]);i++){
    double numWeight = weight[i];
    return numWeight;
  }
  return numWeight*sqrt(stanDev/numData)/numData;
}
