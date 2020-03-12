/*
 *  TDISMC.cpp
 *  
 *  Created by K.Park on 08/10/2016
 *  Copyright 2016. All rights reserved.
 * 
 * This code is for ep or eD : by change ABeam in "TDISMC_EIC.h" 
 */

// General definition in this header
#include "TDISMC_EIC.h"
#include "tim_hobbs/timhobbs.h" //! no need to call because parametrization by fit
#include <time.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH3F.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <TROOT.h>

// init DIS cteq pdf
void initcteqpdf();

int getIndexOfLargestElement(Double_t arr[], int FromIndex, int ToIndex);

// Define form factor function for proton and pion : f2p and f2pi to get TDIS cross-section
// all belows are beased on "neutral pion": DISSOCIATION (dis=1, p -> p + pi0)
double f2p(double x);
double f2pi(double p, double x, double th);
double f2picov(double p, double x, double th);
double f2pitexp(double p, double x, double th);
double f2pipaulvillas(double p, double x, double th);

// Define the cross-section
double F2N(double x, double q2, int inucl);
// new definition for collider frame to call CTEQ function
double cdissigma(double x, double y, double q2, double nu, double ep, int inucl );
double cdissigma_n(double x, double y_D, double q2, double nu, double eprime);
double cdissigma_p(double x, double y_D, double q2, double nu, double eprime);


// Definition of Tim Hobbs functions
double GRV_xVpi(double x,double q2);
double GRV_xSpi(double x,double q2);

// Define function calls for beam smearing
double sigma_th(double pInc, double mInc, double NormEmit, double betaSt);

int mainx(double xMin,double xMax, double Q2Min,double Q2Max, double rnum, const int nevts, const double pbeam, const double kbeam){

  NEvts = nevts;
  PBeam = pbeam;
  kBeam = kbeam;
  
  int e_particle_charge,d_particle_charge,sp_particle_charge,pi_particle_charge,pr_particle_charge;
  e_particle_charge=-1;
  d_particle_charge=1;
  sp_particle_charge=1;
  pr_particle_charge=1;

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  initcteqpdf();
  
  TRandom3 ran3;// Random number generator [45ns/call]
  ran3.SetSeed(rnum);// Sets random generator seed which is used initialize the random number generator
  
  Bool_t iran=kTRUE;      // TRUE==include incident beam emittance
  //Bool_t iran=kFALSE;      // Falut NO incident beam emittance
  
  Bool_t iProton = kTRUE; // TRUE == Spectator Proton
  //Bool_t iProton = kFALSE; // FALSE == Spectator Neutron

  double MSpectator,mPion,MNucleon1;
  int NumPtls,inucl,pi_particle_id;
  double W2;
  
  if(ABeam>1.){
    NumPtls =5; // number of particle in final state (e, D, e', p, ps) // no need to detect pion
    MNucleon1 = MNeutron; 
    mPion = mPim;  // pi-
    pi_particle_id = -211;
    pi_particle_charge=-1;
    inucl = 2; // for neutron target in deuteron
  }
  else{
    NumPtls =5; // Ion beam is proton 
    MNucleon1 = MProton; 
    mPion = mPi0; // pi0
    pi_particle_id = 111;
    pi_particle_charge=0;
    inucl = 1; 
  }
  
  int charge;

  double pprx_inc,ppry_inc,pprz_inc,EprE_inc ;
  
  // STEG spectator proton (always) information
  double_t psx_Lab,psy_Lab,psz_Lab,EsE_Lab ;
  double_t psx_Res,psy_Res,psz_Res,EsE_Res ;
  double_t vsx_Lab,vsy_Lab,vsz_Lab ;
  double_t vsx_Res,vsy_Res,vsz_Res ;
  int sp_particle_id;
  double spmass;
  //    double alpha_r;

  // TDIS scattered electron information
  double pex_Res,pey_Res,pez_Res,EeE_Res;
  double pex_Lab,pey_Lab,pez_Lab,EeE_Lab ;
  double_t vex_Lab,vey_Lab,vez_Lab ;
  double_t vex_Res,vey_Res,vez_Res ;
  int e_particle_id;
  double emass;
  emass =  mElectron;
  e_particle_id = 11;
 
  // TDIS detected pion information
  double ppix_Res,ppiy_Res,ppiz_Res,EpiE_Res;
  double ppix_Lab,ppiy_Lab,ppiz_Lab,EpiE_Lab ;
  double_t vpix_Lab,vpiy_Lab,vpiz_Lab ;
  double_t vpix_Res,vpiy_Res,vpiz_Res ;
  double pimass;
  pimass =  mPion;
 
  // TDIS spectator proton information
  double pprx_Lab,ppry_Lab,pprz_Lab,EprE_Lab ;
  double_t vprx_Lab,vpry_Lab,vprz_Lab ;
  int pr_particle_id;
  double prmass;
  prmass =  MProton;
  pr_particle_id =  2212;

  double xVpiR, xSpiR;
  double yy,rho_PS;
  double F2pi,kTk,pi_int,F2piK;

  // Nuclear codes are given as 10-digit numbers +/-10LZZZAAAI
  // For a (hyper)nucleus consisting of 'n_p' protons, 'n_n' neutrons and
  // 'n_lam' "lam"s, A = n_p + n_n + n_lam gives the total baryon number, 
  // Z = n_p :    the total charge
  // L = n_lam :  the total number of strange quarks.
  // I :          gives the isomer level, 
  // with I = 0 corresponding to the ground state and I > 0 to excitations[1],
  //  where states denoted m, n, p, q translate to I = 1-4. 
  // [1] G. Audi et al., Nucl. Phys. A729, 3 (2003)
  int d_particle_id;
  d_particle_id = 1000010020;
  
  // vertex for the initial particles
  double_t ve0X_Lab,ve0Y_Lab,ve0Z_Lab ;
  double_t vd0X_Lab,vd0Y_Lab,vd0Z_Lab ;

  double weight_tdis;

  // missing particle (X) information
  double pXx_Res,pXy_Res,pXz_Res,EXE_Res ;
  double_t vXx_Res,vXy_Res,vXz_Res ;
  double pXx_Lab,pXy_Lab,pXz_Lab,EXE_Lab ;
  double_t vXx_Lab,vXy_Lab,vXz_Lab ;
  double dummass,zeros;
  int dum_particle_id;
  dum_particle_id = 71;

  if (iProton) {
    MSpectator=MProton;
    charge =1;
    spmass =MSpectator;
    sp_particle_id = 2212;
  }else {
    MSpectator=MNeutron;
    charge =0;
    spmass =MSpectator;
    sp_particle_id = 2112;
  }
  
  char tTitle[80], tName[18], rName[32];
  
  sprintf(tTitle,"p(e,e'np)X Event Generation %3.0f GeV/c x%4.0f GeV/c",kBeam, PBeam);
  sprintf(rName,"../OUTPUTS/pi_p_%.0fon%.0f.root", kBeam, PBeam);

  TFile fRoot(rName,"Recreate", tTitle);
  sprintf(tName,"pi_p");

  TTree *tree = new TTree("Evnts",tTitle);
  
  TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex, kScattered_Vertex;

  // define daughter particle Vertex
  TLorentzVector pSpectator_Vertex,  PX_Vertex, PX_Vertex_Rest;
  TLorentzVector pScatterProton_Vertex;
  TLorentzVector pScatterPion_Vertex;

  typedef struct{
    Double_t s_e, s_q,  Q2, xBj, nu, W, p_RT, tempVar,
      pDrest, x_D, y_D, Yplus,  tSpectator, tPrime,
      TwoPdotk, TwoPdotq, MX2,
      alphaS, pPerpS, pPerpZ, Jacob;
  } Invariants;
  static Invariants invts;

  Double_t Jacob;

  //  ************************************************
  // Store physics variable for TDIS into ROOTs
  //  ************************************************
  double tpi, ypi, y, fpi, xpi; 
  double sigma_dis;
  double sigma_tdis;
  double f2N;   
  double p2_pt, p2_z;
  Double_t qMag, pDotq;
	
  double TDIS_xbj, TDIS_znq, TDIS_Mx2, TDIS_y;

  const Int_t bufsize=32000;

  // invariants
  tree->Branch("invts",&invts,		     "s_e/D:s_q/D:Q2/D:xBj/D:nu/D:W/D:x_D/D:y_D/D:Yplus/D:tSpectator/D:tPrime/D:TwoPdotk/D:TwoPdotq/D:p_RT/D:pDrest/D:tempVar/D:MX2/D:alphaS/D:pPerpS/D:pPerpZ/D");
  // incoming particle electron, deuteron(proton)
  tree->Branch("e_Inc.", &kIncident_Vertex,bufsize,1);
  tree->Branch("P_Inc.",  &PIncident_Vertex, bufsize, 1);
  // outgoing particles (electron, virtual photon, proton1, proton2, pion)
  tree->Branch("e_Scat.", &kScattered_Vertex, bufsize, 1);
  tree->Branch("q_Vir.", &qVirtual_Vertex, bufsize, 1);
  tree->Branch("p1_Sp.", &pScatterProton_Vertex, bufsize, 1);
  tree->Branch("p2_Sp.", &pSpectator_Vertex, bufsize, 1);
  tree->Branch("pi.", &pScatterPion_Vertex, bufsize, 1);
  //	tree->Branch("psf",&psf,"psf/D");
  //tree->Branch("Jacob",&Jacob,"Jacob/D");

  // Store physics varaibles into ROOTs
  tree->Branch("xpi", &xpi, "xpi/D");
  tree->Branch("ypi", &ypi, "ypi/D");
  tree->Branch("y", &y, "y/D");
  tree->Branch("tpi", &tpi, "tpi/D");
  tree->Branch("fpi", &fpi, "fpi/D");
	
  tree->Branch("f2N", &f2N, "f2N/D");
  tree->Branch("p2_pt", &p2_pt, "p2_pt/D");
  tree->Branch("p2_z", &p2_z, "p2_z/D");
	
  tree->Branch("sigma_dis", &sigma_dis, "sigma_dis/D");
  tree->Branch("sigma_tdis", &sigma_tdis, "sigma_tdis/D");

  // Incident proton
  tree->Branch("pprx_inc", &pprx_inc, "pprx_inc/D");
  tree->Branch("ppry_inc", &ppry_inc, "ppry_inc/D");
  tree->Branch("pprz_inc", &pprz_inc, "pprz_inc/D");
  tree->Branch("EprE_inc", &EprE_inc, "EprE_inc/D");

  // TDIS scattered electron information
  tree->Branch("psx_Lab", &psx_Lab, "psx_Lab/D");
  tree->Branch("psy_Lab", &psy_Lab, "psy_Lab/D");
  tree->Branch("psz_Lab", &psz_Lab, "psz_Lab/D");
  tree->Branch("EsE_Lab", &EsE_Lab, "EsE_Lab/D");
  
  // TDIS scattered electron information
  tree->Branch("pex_Lab", &pex_Lab, "pex_Lab/D");
  tree->Branch("pey_Lab", &pey_Lab, "pey_Lab/D");
  tree->Branch("pez_Lab", &pez_Lab, "pez_Lab/D");
  tree->Branch("EeE_Lab", &EeE_Lab, "EeE_Lab/D");

  // TDIS detected pion information
  tree->Branch("ppix_Lab", &ppix_Lab, "ppix_Lab/D");
  tree->Branch("ppiy_Lab", &ppiy_Lab, "ppiy_Lab/D");
  tree->Branch("ppiz_Lab", &ppiz_Lab, "ppiz_Lab/D");
  tree->Branch("EpiE_Lab", &EpiE_Lab, "EpiE_Lab/D");

  // TDIS spectator proton information
  tree->Branch("pprx_Lab", &pprx_Lab, "pprx_Lab/D");
  tree->Branch("ppry_Lab", &ppry_Lab, "ppry_Lab/D");
  tree->Branch("pprz_Lab", &pprz_Lab, "pprz_Lab/D");
  tree->Branch("EprE_Lab", &EprE_Lab, "EprE_Lab/D");

  // missing particle (X) information
  tree->Branch("pXx_Lab", &pXx_Lab, "pXx_Lab/D");
  tree->Branch("pXy_Lab", &pXy_Lab, "pXy_Lab/D");
  tree->Branch("pXz_Lab", &pXz_Lab, "pXz_Lab/D");
  tree->Branch("EXE_Lab", &EXE_Lab, "EXE_Lab/D");
	
  tree->Branch("TDIS_xbj", &TDIS_xbj, "TDIS_xbj/D");
  tree->Branch("TDIS_znq", &TDIS_znq, "TDIS_znq/D");
  tree->Branch("TDIS_Mx2", &TDIS_Mx2, "TDIS_Mx2/D");
  tree->Branch("TDIS_y", &TDIS_y, "TDIS_y/D");
	
  double sig_eThx=0.0, sig_eThy=0.0, sig_iThx=0.0, sig_iThy=0.0;
  //	double sig_k=0.0, sig_P=0.0;
  double MIon;
  // Global Phase Space Factor;
  //s_e - M_ion^2 - m_electron^2
  double uu, uv, uw, ux, uy, norm;
  
  if (ZBeam==1.0&&ABeam==2.0) {
    MIon = MDeut;
  } else {
    MIon = MProton;
  }
  
  double PIon = PBeam*ZBeam;
  double s_e0 = 2.*kBeam*(sqrt(PIon*PIon+MIon*MIon)+PIon*cos(CrossingTheta))+ MIon*MIon + mElectron*mElectron;
  
  printf("Your kinematics: [xBj_min:xBj_max] = [%9.6f:%9.6f] \n", xMin,xMax);
  printf("Your kinematics: [Q2_min:Q2_max] = [%9.6f:%9.6f] \n", Q2Min,Q2Max);
  printf("Incident Ion Mass %9.5f GeV \n", MIon);
  printf("Incident Electron, Ion Momenta: %8.4f, %8.2f GeV/c | s_0 = %10.4f GeV^2 \n", kBeam, PIon, s_e0);   
	
  if (iran) {
    sig_eThx = sigma_th(kBeam, mElectron, eEpsNX, eBetaStarX);
    sig_eThy = sigma_th(kBeam, mElectron, eEpsNY, eBetaStarY);
    sig_iThx = sigma_th(PBeam, MIon, iEpsNX, iBetaStarX);
    sig_iThy = sigma_th(PBeam, MIon, iEpsNY, iBetaStarY);
  }
  
  double kBeamMC, kBeamMCx, kBeamMCy, kBeamMCz;
  double PBeamMC, PBeamMCx, PBeamMCy, PBeamMCz;
  double EScatRest, kScatRest, csTheRest, PhiScatRest, csPhiRest;
  double EScatRest0, kScatRest0, csTheRest0, PhiScatRest0;

  tree->Branch("EScatRest", &EScatRest, "EScatRest/D");
  tree->Branch("kScatRest", &kScatRest, "kScatRest/D");
  tree->Branch("PhiScatRest", &PhiScatRest, "PhiScatRest/D");
  tree->Branch("csPhiRest", &csPhiRest, "csPhiRest/D");
  tree->Branch("csTheRest", &csTheRest, "csTheRest/D");
	
  TVector3       UnitXLab, UnitYLab, UnitZLab;
  
  UnitXLab.SetXYZ(1.0,0.0,0.0);
  UnitYLab.SetXYZ(0.0,1.0,0.0);
  UnitZLab.SetXYZ(0.0,0.0,1.0);
  
  TVector3       UnitXRest, UnitYRest, UnitZRest;
  TVector3       UnitXqCM,  UnitYqCM,  UnitZqCM;
  TVector3       BoostCM, BoostRest, BoostRest0;
  TVector3       kScat3vec, pS_3vec;
  TVector3       kScat3vec0, pS_3vec0;
  
  TLorentzVector kIncident_Rest, PIncident_Rest, qVirtual_Rest, p_DT;
  TLorentzVector kIncident_Rest0, PIncident_Rest0, qVirtual_Rest0, p_DT0, p_ST0, p_ST;
  TLorentzVector kScattered_Rest, pSpectator_Rest, PX_Vector_Rest;
  TLorentzVector kScattered_Rest0, pSpectator_Rest0, PX_Vector_Rest0;
  TLorentzVector kIncident_0, PIncident_0; // unsmeared
  
  kIncident_0.SetXYZM(0.0,0.0,-kBeam,mElectron);
  PIncident_0.SetXYZM(PIon*sin(CrossingTheta)*cos(CrossingPhi),PIon*sin(CrossingTheta)*sin(CrossingPhi),PIon*cos(CrossingTheta),MIon);
  
  // TTree *itree = new TTree("Init",tTitle);
  // itree->Branch("e0.", &kIncident_0,bufsize,1);
  // itree->Branch("P0.",  &PIncident_0, bufsize, 1);
  
  typedef struct {
    Int_t nEvts;
    Float_t PhSpFct;
  } MonteCarlo;
  static	MonteCarlo mc;
  mc.nEvts = NEvts;
  mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow((pSMax),3)/3.;
  // itree->Branch("MC", &mc, "nEvts/I:PhSpFct/F");

  double pS_rest, csThRecoil, phiRecoil;

  //name of output file : = "TDIS_lund.dat";
  ofstream OUT (Form("../OUTPUTS/pi_p_lund_%.0fon%.0f.dat", kBeam, PBeam), ios::app);

  // **********************************************************************************
  // define TDIS pSpectator with fermi momentum from  data file "moment_ld2b.dat" from G4SBS 
  // **********************************************************************************
  // This is initial neutron momentum distribution with fermi motion (from initial deuteron beam)
  TVector3 V3_p1,V3n_p1;
  double V3Px_p1, V3Py_p1, V3Pz_p1;
  double Px_p1, Py_p1, Pz_p1, P_p1, E_p1, theta_p1, phi_p1;

  // Re-defined the LorentzVector after taking into count fermi momentum : REST frame
  TLorentzVector pSpectator_RestNew;

  double phi_p2;
  
  // TDIS-SBS xbj means xbj_D (xbj of nucleon should be xbj/2 ? Asking Dasuni for xbj definition
  double Px_p2,Py_p2,Pz_p2,P_p2,theta_p2,E_p2;

  // Get LorentzVector for the pScattered Proton for TDIS in rest frame
  TLorentzVector pScatterProton_Rest;

  double ppr_Lab;

  // *****************************************************************************
  // Generate missing Hadron(Pion) for TDIS in rest frame
  // initial neutron 4momentum(-P1) - TDIS scattered proton 4momentum(P2)
  // *****************************************************************************
  double E_pi,Px_pi,Py_pi,Pz_pi;

  TLorentzVector pScatterPion_Rest;

  TVector3 pScatterPion_V3;

  double P_pi;

  double fpifac;

  double minE = 0.0;
  
  //for progress bar
  double progress=0.0;

  Int_t MEvts=1;
  
  // Starting of event creating
  for (int iEvt=0; iEvt<=NEvts; iEvt++) {

    if(iEvt%1000==0) {	    
      int barWidth = 70;
      progress = ((double)iEvt/(double)NEvts);	    
      // cout<<iEvt<<"/"<<NEvts<<endl;
      // cout << progress << endl;
      std::cout << "[";
      double pos = barWidth * progress;
      for (double i = 0.; i < barWidth; ++i) {
	if (i < pos) std::cout << "=";
	else if (i == pos) std::cout << ">";
	else std::cout << " ";
      }
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();	    
    }	 
    
    Jacob = 1.0;
	  
    if (iran) {  // taking into account initial beam smearing
      kBeamMC = kBeam*ran3.Gaus(1.0,eDkOverk);
      kBeamMCx= kBeamMC*ran3.Gaus(0.0,sig_eThx);
      kBeamMCy= kBeamMC*ran3.Gaus(0.0,sig_eThy);
      kBeamMCz=-kBeamMC;  //  Angles are really Tangents
      PBeamMC = PIon*ran3.Gaus(1.0,iDPoverP);
      PBeamMCx= PBeamMC*ran3.Gaus(0.0,sig_iThx); // IP1 configuration
      PBeamMCy= PBeamMC*ran3.Gaus(0.0,sig_iThy);
      PBeamMCz= PBeamMC;	    
    } else { // no beam smearing
      kBeamMC = kBeam;
      kBeamMCx= 0.0;
      kBeamMCy= 0.0;
      kBeamMCz= -kBeam;
      PBeamMC = PIon;
      PBeamMCx= 0.0;
      PBeamMCy= 0.0;
      PBeamMCz= PBeamMC;
    }
    
    //  definition of beam particle vertex
    ve0X_Lab = 0;
    ve0Y_Lab = 0;
    ve0Z_Lab = 0;
    vd0X_Lab = 0;
    vd0Y_Lab = 0;
    vd0Z_Lab = 0;	 	
		
    // these are with beam smearing
    kIncident_Vertex.SetXYZM(kBeamMCx, kBeamMCy, kBeamMCz, mElectron);// Set 4-momentum of incident electron beam
    PIncident_Vertex.SetXYZM(PBeamMCx, PBeamMCy, PBeamMCz, MIon);// Set 4-momentum of incident ion beam
    
    //  Crossing angle JLEIC with smearing
    PIncident_Vertex.RotateY(CrossingTheta);
    PIncident_Vertex.RotateZ(CrossingPhi);
		
    // Calculation of invariant variable from two beam particles (s_e, s_q) with smearing
    invts.TwoPdotk = 2.*PIncident_Vertex.Dot(kIncident_Vertex);
    invts.s_e = MIon*MIon + mElectron*mElectron + invts.TwoPdotk;
    invts.TwoPdotq = 2.*PIncident_Vertex.Dot(qVirtual_Vertex); 
    invts.s_q = MIon*MIon + invts.TwoPdotq;   
	  
    // Boosting vectors from electron+Ion CM frame
    // Boosting vector from Ion rest frame
    BoostRest = PIncident_Vertex.BoostVector(); 		  
	  
    PIncident_Rest = PIncident_Vertex; 

    PIncident_Rest.Boost(-BoostRest); // should result in (0.,0.,0.,MIon)
    kIncident_Rest = kIncident_Vertex;
    kIncident_Rest.Boost(-BoostRest);

    pprx_inc = PIncident_Vertex.X();
    ppry_inc = PIncident_Vertex.Y();
    pprz_inc = PIncident_Vertex.Z();
    
    EprE_inc = PIncident_Vertex.E();
    // EprE_inc = sqrt(pprx_inc*pprx_inc+ppry_inc*ppry_inc+pprz_inc*pprz_inc+MProton*MProton);

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    /*
      cout << " ---------------------------------------------------------------------------------------- "  << endl;
      cout << " ---------------------------------------------------------------------------------------- "  << endl;
      cout << "(A) kIncident_Vertex (Lab)= " << kIncident_Vertex.X() << ", " << kIncident_Vertex.Y() << ", " << kIncident_Vertex.Z() << endl;
      cout << "(B) kIncident_Rest (rest)= " << kIncident_Rest.X() << ", " << kIncident_Rest.Y() << ", " << kIncident_Rest.Z() << endl;
    //*/
		
    // We generated xBj and Q2 binning...
    // Q is not calculated by Ei- Ef
    // Ef is calculated by given Q2 from random number
	
    // Generate Q2, ln(xBj) uniformly
    //  psf *= (Q2Max-Q2Min)*(log(xMax)-log(xMin));
    
    uu   = ran3.Uniform();
    invts.Q2   = Q2Max*uu + Q2Min*(1.-uu);
    uv   = ran3.Uniform();
    invts.xBj  = pow(xMin,1.-uv)*pow(xMax,uv);
    invts.x_D  = invts.xBj*(MProton/MIon);
    invts.y_D  = invts.Q2/(invts.x_D*invts.TwoPdotk);
    invts.Yplus = 1 + ((1-invts.y_D)*(1-invts.y_D));

    y = invts.y_D;

    if (invts.y_D>=(1.0-2.*mElectron*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
    if (invts.y_D>=1./(1.+invts.x_D*MIon*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
	  
    Jacob *=invts.xBj;   // Jacobian  dx/dlnx
    //		psf *=2.*pi;
    // scattered electron, phi angle in rest frame 
    PhiScatRest = pi*(2.*ran3.Uniform()-1.0);
    csPhiRest = cos(PhiScatRest);
    //  mElectron-->0 approximation
    EScatRest = invts.TwoPdotk*(1.-invts.y_D)/(2.*MIon);
    
    if (EScatRest<mElectron) {
      // should never happen
      printf("illegal Rest frame scattered electron energy =%10.6f \n",EScatRest);
      continue;
    }
    
    kScatRest = sqrt(EScatRest*EScatRest-mElectron*mElectron);// scattered electron momentum in rest frame
    csTheRest = (2.*EScatRest*kIncident_Rest.E() - invts.Q2 - 2.*mElectron*mElectron) /(2.*kScatRest*kIncident_Rest.P());   // scattered electron \cosine\theta in rest frame	                                                 	  
		
    if (csTheRest*csTheRest>1.0) {
      // should never happen
      printf("illegal Rest frame cos(the) = %6.2f \n", csTheRest);
      printf(" (k_Rest, k'_Rest, Q2, xBj) = (%8.3f,  %8.4f, %6.2f, %5.3f) \n",kIncident_Rest.E(),kScatRest,invts.Q2,invts.xBj);
      printf(" (2k.P, invts.s_e, y_D, x_D) = (%6.2f, %6.2f,%10.6f, %8.4f) \n",invts.TwoPdotk, invts.s_e, invts.y_D, invts.x_D);
      continue;
    }
	  
    // Definition of unit vector in rest frame (norm vector)
    UnitZRest  = -kIncident_Rest.Vect();
    norm       = 1./UnitZRest.Mag();
    UnitZRest *= norm;
    UnitYRest  = UnitZRest.Cross(UnitXLab);
    norm       = 1./UnitYRest.Mag();
    UnitYRest *= norm;
    UnitXRest  = UnitYRest.Cross(UnitZRest);
		
    // 3-vector of scattered electron in rest frame with smearing
    kScat3vec = kScatRest*(sin(acos(csTheRest))*(cos(PhiScatRest)*UnitXRest+sin(PhiScatRest)*UnitYRest)-csTheRest*UnitZRest);
    kScattered_Rest.SetVectM(kScat3vec, mElectron);
    qVirtual_Rest = kIncident_Rest-kScattered_Rest;
    kScattered_Vertex = kScattered_Rest;		

    // *** Scattered electron in REST frame
    // For debugging purpose: checked
    /*
      cout << kScat3vec(0) << "  " << kScat3vec(1) << "  " << kScat3vec(2)  << " sqrt(x*x+y*y+z*z) =  " << kScatRest << endl;
    */
	  
    // Deuteron momentum in rest frame with and without smearing
    invts.pDrest = sqrt(PIncident_Rest(0)*PIncident_Rest(0)+PIncident_Rest(1)*PIncident_Rest(1)+PIncident_Rest(2)*PIncident_Rest(2));

    // ***
    // for the debugging purpose: OK, CHECKED
    /*
      cout << "electron rest frame =" << kScat3vec.X() << ", " << kScat3vec.Y() << ", " << kScat3vec.Z() << endl;
    */
	  
    // Back to Lab frame
    kScattered_Vertex.Boost(BoostRest);
    qVirtual_Vertex = kIncident_Vertex-kScattered_Vertex;
	  
    // Hadronic Unit vectors relative to q
    UnitZqCM = -qVirtual_Rest.Vect();
    norm     = 1./UnitZqCM.Mag();
    UnitZqCM*= norm;
    UnitYqCM = -kScat3vec.Cross(kIncident_Rest.Vect());
    norm     = 1./UnitYqCM.Mag();
    UnitYqCM*=norm;
    UnitXqCM = UnitYqCM.Cross(UnitZqCM);

    qMag = qVirtual_Rest.P();
    pDotq = BoostRest(0)*qVirtual_Rest(0)+BoostRest(1)*qVirtual_Rest(1)+BoostRest(2)*qVirtual_Rest(2);
    p_DT = PIncident_Vertex - (pDotq/qMag/qMag)*qVirtual_Rest;
    p_ST = pSpectator_Vertex -(pDotq/qMag/qMag)*qVirtual_Rest;		
    
    V3Px_p1    = -PIncident_Rest.X()/ABeam; // simply opposite of fermi momentum in rest frame (same as TDIS SBS)
    V3Py_p1    = -PIncident_Rest.Y()/ABeam; // n,p momentum carrying : half of mementum of deuteron
    V3Pz_p1    = -PIncident_Rest.Z()/ABeam;
    V3_p1.SetXYZ(V3Px_p1,V3Py_p1,V3Pz_p1);
    V3n_p1    = V3_p1 + V3_p1.Unit();  // adding up fermi momentum
    
    // ***
    // For debugging purpose: OK
    /*
      cout << " ---------------------------------------------------------------------------------------- "  << endl;
      cout << "(0) testing print init neutron momentum in rest ::: " << V3n_p1.X() << ",  " << V3n_p1.Y() <<  ",  " << V3n_p1.Z() << endl;
    */
    
    if(ABeam>1.){
      Px_p1    = V3n_p1.X();
      Py_p1    = V3n_p1.Y();
      Pz_p1    = V3n_p1.Z();
    }
    else{  // free proton no fermi momentum in REST frame
      Px_p1    = 0.;
      Py_p1    = 0.;
      Pz_p1    = 0.;
    }
		
    P_p1     = V3n_p1.Mag();
    E_p1     = sqrt (pow(P_p1,2) + pow(MNeutron,2));
    theta_p1 = acos( Pz_p1/P_p1);

    phi_p1 = 0.0;
    
    if (Py_p1 > 0.0){
      phi_p1 = acos( Px_p1/sqrt(Px_p1*Px_p1 + Py_p1*Py_p1));}
    else {
      phi_p1 = 2*PI - acos(Px_p1/sqrt(Px_p1*Px_p1 + Py_p1*Py_p1));
    }
    
    // *** This proton is TDIS initial neutron with fermi momentum 
    // for the debugging purpose
    /*
      cout << "(1) Beam initial neutron from D in RestFrame =" << V3Px_p1 << ", " << V3Py_p1 << ", " << V3Pz_p1 << ", "  << endl;
      cout << "(2) Beam initial neutron from D with fermi momentum in RestFrame =" << Px_p1 << ", " << Py_p1 << ", " << Pz_p1 << ", "<< "should be same as (0)" << endl;
      cout << "(***) invts.TwoPdotk= " <<  invts.TwoPdotk << ", invts.y_D= " << invts.y_D << ", EScatRest= " << EScatRest << endl;
    */
    
    pSpectator_RestNew.SetXYZM(Px_p1,Py_p1,Pz_p1,MNucleon1);
		
    // *****************************************************************************
    // Generate pScattered Proton for TDIS with Spherically Symmetric in rest frame (this is P2: scattered proton)
    // neutron -> proton + pi
    // *****************************************************************************
    // TDIS in SBS software, they define the speactor(proton) with pt[0.:0.5GeV], z[0.:1.0] with plat random
    // in TDIS kinematics, the general Spectator particle (proton) is nothing to do with following calculation.
    //
    // neutron config = proton + pion(-):  this proton called additional speactator in TDIS concept
    // randomize in pt and z
    p2_pt = gRandom->Uniform(0.005*PBeam); // .5% of incoming ion beam momentum, this is the limit of the transverse momentum of  recoil particle...
    p2_z = gRandom->Uniform(1.);

    //  definition are moved at the beginning of code
    //		double TDIS_xbj, TDIS_znq,TDIS_Mx2,TDIS_y;
    TDIS_xbj = invts.Q2/(2*pSpectator_RestNew.Dot(qVirtual_Rest));
    TDIS_znq = p2_z*pSpectator_RestNew.Dot(qVirtual_Rest);
    TDIS_Mx2 = (qVirtual_Rest + pSpectator_RestNew).Mag2();
    TDIS_y   = (pSpectator_RestNew.Dot(qVirtual_Rest))/(pSpectator_RestNew.Dot(kIncident_Rest));
		
    //***
    // For debugging purpose:
    /*
      cout << "(3) kinematics: TDIS.xBj = " << TDIS_xbj <<  ", invts.xBj= " << invts.xBj <<  ", invts.Q2= " << invts.Q2<< endl; 
      cout << "(3*) generating: p2_pt " << p2_pt <<  ", p2_z= " << p2_z <<  ", phi_p2= " << phi_p2 << endl; 
    */

    phi_p2 = gRandom->Uniform(360.0*D2R);

    Pz_p2 = ( -1.0*TDIS_znq*qVirtual_Rest.Z() +sqrt( TDIS_znq*qVirtual_Rest.Z()*TDIS_znq*qVirtual_Rest.Z() +invts.Q2*(qVirtual_Rest.E()*qVirtual_Rest.E())*(MProton*MProton + p2_pt*p2_pt)-invts.Q2*TDIS_znq*TDIS_znq) )/invts.Q2/1000.;  // unit match with "p2_pt" (GeV)
	  
    P_p2  = sqrt (Pz_p2*Pz_p2 + p2_pt*p2_pt);
    theta_p2 = acos (Pz_p2/P_p2);
    E_p2 = sqrt ( P_p2*P_p2 + MProton*MProton);
	  
    Px_p2 = p2_pt*cos(phi_p2);
    Py_p2 = p2_pt*sin(phi_p2);

    // for dubugging purpose:
    // Generated additional spectator proton from TDIS reaction
    /*
      cout << "\n" << "(4**) Spectator(T), Pz_p2= " << Pz_p2 << ", P_p2 =" << P_p2  << endl;
    */
	  
    // for the debugging purpose: SEEMS NOT CRAZY NUMBER....
    //cout << "TDIS missing mass =" << TMath::Sqrt(TDIS_Mx2)  << endl;
	  
    pScatterProton_Rest.SetXYZM(P_p2*sin(theta_p2)*cos(phi_p2), P_p2*sin(theta_p2)*sin(phi_p2), P_p2*cos(theta_p2),MProton);
    
    // *****************************************************************************

    // for dubugging purpose:
    // Generated additional spectator proton from TDIS reaction
    /*
      cout << "\n" << "(4) Spectator(T), px= " << pScatterProton_Rest.X() << ", py= " << pScatterProton_Rest.Y()<<
      ", pz= " << pScatterProton_Rest.Z() << ", it's random numbers (pt,z)"<< endl;	    
      cout  << "(5) Active(N), px= " << pSpectator_RestNew.X() << ", py= " << pSpectator_RestNew.Y()<<
      ", pz= " << pSpectator_RestNew.Z() << ",  Should be same as (2)"<< endl;	    
      cout  << "(6) qVirtual REST, px= " << qVirtual_Rest.X() << ", py= " << qVirtual_Rest.Y()<<
      ", pz= " << qVirtual_Rest.Z() << endl;
    */

    // Back to Lab frame
    pScatterProton_Vertex = pScatterProton_Rest;
    pScatterProton_Vertex.Boost(BoostRest);
    pprx_Lab = pScatterProton_Vertex.X();
    ppry_Lab = pScatterProton_Vertex.Y();
    pprz_Lab = pScatterProton_Vertex.Z();
    EprE_Lab = sqrt(pprx_Lab*pprx_Lab+ppry_Lab*ppry_Lab+pprz_Lab*pprz_Lab+MProton*MProton);
    vprx_Lab = 0.0;
    vpry_Lab = 0.0;
    vprz_Lab = 0.0;

    ppr_Lab = sqrt(pprx_Lab*pprx_Lab+ppry_Lab*ppry_Lab+pprz_Lab*pprz_Lab);
	  
    // For debugging purpose
    /*
      if (ppr_Lab > (PBeam/ABeam) ) { 
      // Unphysical kinematics  // if TDIS spectator momentum larger than 50% ion momentum
      printf("impossible of TDIS spectator momentum= %6.2f \n", ppr_Lab);
      continue;
      }
    */
    
    E_pi  = pSpectator_RestNew.E() - pScatterProton_Rest.E();
    Px_pi = pSpectator_RestNew.X() - pScatterProton_Rest.X(); 
    Py_pi = pSpectator_RestNew.Y() - pScatterProton_Rest.Y();
    Pz_pi = pSpectator_RestNew.Z() - pScatterProton_Rest.Z();
		
    pScatterPion_Rest.SetXYZM(Px_pi,Py_pi,Pz_pi,mPion);
    
    pScatterPion_V3.SetXYZ(Px_pi,Py_pi,Pz_pi);
		
    // Back to Lab frame
    pScatterPion_Vertex = pScatterPion_Rest;
    pScatterPion_Vertex.Boost(BoostRest);
    ppix_Lab = pScatterPion_Vertex.X();
    ppiy_Lab = pScatterPion_Vertex.Y();
    ppiz_Lab = pScatterPion_Vertex.Z();
    EpiE_Lab = sqrt(ppix_Lab*ppix_Lab+ppiy_Lab*ppiy_Lab+ppiz_Lab*ppiz_Lab+mPion*mPion);

    // for debuggin purpose
    /*
      cout  << "(7L) pion Vertex, px= " << ppix_Lab << ", py= " << ppiy_Lab <<
      ", pz= " << ppiz_Lab << ", Epi= " << EpiE_Lab << endl;
    */
    
    vpix_Lab = 0.0;
    vpiy_Lab = 0.0;
    vpiz_Lab = 0.0;
		
    // Definition of pion variable
    //double tpi, ypi, fpi, xpi; 
    xpi = TDIS_xbj/(1 - p2_z);
    tpi = (E_pi*E_pi) - (pScatterPion_V3.Mag()*pScatterPion_V3.Mag());
    ypi = pScatterPion_Rest.Dot(qVirtual_Rest)/(pScatterPion_Rest.Dot(kIncident_Rest));
	  
    P_pi = pScatterPion_V3.Mag();
	  
    // *****
    // for debugging purpose
    /*
      cout  << "(7R) pion REST, px= " << pScatterPion_Rest.X() << ", py= " << pScatterPion_Rest.Y()<<
      ", pz= " << pScatterPion_Rest.Z() << ", E= " << pScatterPion_Rest.E() << ", M= " << pScatterPion_Rest.M() << endl;
    */
    
    if(ABeam>1){
      fpifac = 2;
    }
    else{
      fpifac = 1;
    }
	  
    // pion form factor is only defined in xBj = (0.055 : 0.3) --- Removed this constraint in TDISMC_EIC.h
    // if (TDIS_xbj > 0.001 && TDIS_xbj < 1.0){
    if (xpi > 0.001 && xpi < 1.0){ // xpi for ZEUS parameterization

      //  typ = 2 ! s-dependent expon FORM FACTOR (L=1.56,dis=1,FLAG=0)
      // fpi = fpifac*f2pi(P_pi, TDIS_xbj, theta_p2/D2R);
	    
      //  typ = 3 ! COV DIP FORM FACTOR
      // fpi = fpifac*f2picov(P_pi, TDIS_xbj, theta_p2/D2R);

      //  typ = 5 ! t-dependent expon FORM FACTOR (L=0.85,dis=1,FLAG=0)
      // fpi = fpifac*f2pitexp(P_pi, TDIS_xbj,theta_p2/D2R);
	    
      //  typ = 6 ! Paul-Villas FORM FACTOR (L=0.27,dis=1,FLAG=0)
      // fpi = fpifac*f2pipaulvillas(P_pi, TDIS_xbj, theta_p2/D2R);

      //  ZEUS parameterization with F2N (proton SF)
      // fpi = fpifac*(0.361)*F2N(xpi, invts.Q2, inucl);
      fpi = fpifac*f2piZEUS(TDIS_xbj);

      /*
	L = renormalization cut-off parameter... ex) typ = 2 ! 1.56GeV
	dis = charge state configuration (0=charge exchange, 1=neutral exchange)
	FLAG = combination of splitting function/PDF (0:\pi, 1:\rho, 2:\pi\Delta)
      */

      /*
      // ***************************************************************************************
      // ************************ TESTing with FF of pion exchange model ***********************
      // ***************************************************************************************
	    
      // Testing of "y" variable: assumption: same kinematic variables:  x=TDIS_xbj, k=P_pi, th=theta_p2
      yy = (P_pi/MProton)*cos(theta_p2)+(1/MProton)*(MProton-sqrt(MProton*MProton+P_pi*P_pi));
      // phase space factor (Jacobian)  !
      rho_PS = (2./MProton)*P_pi*P_pi*sin(theta_p2)*(1.0-sin(phi_p2)*cos(theta_p2)); 		  
	    
      // TESTing pion PDF from GRV (inputs: x, Q2) (outputs: xVpi,  xSpi)  *********************
      // by definition (theory paper) x/y = xpi
      xVpiR = GRV_xVpi(xpi,invts.Q2);
      xSpiR = GRV_xSpi(xpi,invts.Q2);
	    		  
      F2pi = (5./9.) *(xVpiR + 2.0* xSpiR);
      kTk = P_pi * sin(theta_p2);
      pi_int = fypiN(yy,kTk,L,typ,dis)*F2pi;

      // |k| integral loop 
      F2piK = (2./3.)*(2./3.)*pi_int;
      if(fpi>0){
      cout << "***************       ******************      *********************" << endl;
      //cout << "test outout GRV xVpi = " << xVpiR << " , at xbj= "<< TDIS_xbj <<", Q2= " << invts.Q2 << endl;
      //cout << "test outout GRV xSpi = " << xSpiR << " , at xbj= "<< TDIS_xbj <<", Q2= " << invts.Q2 << endl;
	    
      }
      // ***************************************************************************************
      // ************************ TESTing with FF of pion exchange model ***********************
      // ***************************************************************************************
      */
		    
    }
    else{
      fpi = 0.0;
    }
		
    // *****************************************************************************

    //fpi = 2*f2pi(P_pi, invts.xBj, theta_p2/D2R);
	  
    // ***
    // For debugging purpose
    /*
      cout << "(8) " <<"fpi= " << fpi <<  ", 2nd proton Pz=" << p2_pt  << ", " << Pz_p2  << ", P_pi= "
      << P_pi << ", xbj= "<< TDIS_xbj << ", theta_p2(deg)= "<< theta_p2 << endl; 
    */

    invts.nu = invts.Q2 / (2 * MIon * invts.xBj) ;
    invts.W = MIon*MIon + invts.Q2/invts.xBj -invts.Q2;
    W2 = invts.W*invts.W;
	      		
    //  **********************************************************
    //
    //    Cross-section assignment
    //
    //  **********************************************************		
    // From the DIS event generator

    if (TDIS_xbj > 0.001 && TDIS_xbj < 1.0){
      // *** for the neutron/proton target
      // **** dissigma(E0, Theta_e, E') : where these are in Lab frame ??? should I apply electron lab ?
      //sigma_dis    = dissigma_n(kIncident_Rest.E(),
      //                          acos( kScattered_Rest.Z()/kScat3vec.Mag()),
      //                          kScattered_Rest.E()) ;
      // *** dissigma_n is also defined in G4SBSDIS.hh
      // f2N = 1.0;
      f2N   = F2N(invts.xBj, invts.Q2, inucl);  // F2N is define at G4SBSDIS.hh
      // sigma_dis = cdissigma(invts.xBj,invts.y_D,invts.Q2,invts.nu,kScattered_Rest.E(),1);
      // *** for the JLEIC collider : it needs an input as (x, y_D, q2, nu, eprime_rest)  instead of ebeam, theta_e, eprime 
      if(ABeam>1.){
	sigma_dis    = cdissigma_n(TDIS_xbj,invts.y_D,invts.Q2,invts.nu,kScattered_Rest.E());
	// sigma_dis    = cdissigma_n(TDIS_xbj,TDIS_y,invts.Q2,invts.nu,kScattered_Rest.E());
      }
      else{
	sigma_dis    = cdissigma_p(TDIS_xbj,invts.y_D,invts.Q2,invts.nu,kScattered_Rest.E());
	// sigma_dis    = cdissigma_p(TDIS_xbj,TDIS_y,invts.Q2,invts.nu,kScattered_Rest.E());
      }
      sigma_tdis   = sigma_dis * (fpi/f2N); 
    // }else{
    //   sigma_dis = 0.0;
    //   sigma_tdis = 0.0;
    //   f2N=0.0;
    }

    /*
      if(fpi>0.){
      cout << " what is f2N= "  << f2N  << ", fpi=" << fpi << ", momentum = " << P_pi << ", theta = " << theta_p2/D2R << endl;
      }
    */
		
    // ***
    // For debugging purpose
    /*
      cout << "(10) " <<"sigma_dis=  " << sigma_dis  <<  ", sigma_tdis= " << sigma_tdis
      << ",  f2N= " <<  f2N << endl;
      cout << "kIncident_Rest.E()= "<< kIncident_Rest.E() << ", acos( kScattered_Rest.Z()/kScattered_Rest.Mag())= "
      <<kScattered_Rest.Z() <<", " << kScat3vec.Mag() << ", kScattered_Rest.E()= " <<kScattered_Rest.E()
      << ", nanobarn= "<< nanobarn<< endl;
      cout << "(11) " <<"invts.nu=  " << qVirtual_Rest.E()/MIon << ", " << invts.nu<< ",  kIncident_Rest.E()= "
      << kIncident_Rest.E() << ",  incoming electron energy= "<< invts.nu*2. + kScattered_Rest.E() << endl; 
    */
	  
    // qVirtual_Rest.E(): the virtual 4- momentum transfer between electron and deuteron in rest frame
    // invts.nu : the virtual 4- momentum transfer between electron and target nucleon in rest frame
	  
    //  invts.nu*2. + kScattered_Rest.E() : initial electron energy in rest frame		

    // If the hadron beam is Ion (deuteron)
	  
    if(ABeam>1.){
	    
      // **********************************************************
      //
      // 
      //                Spectator nucleon for STEG
      //
      //
      // **********************************************************
      // Generate pSpectator Spherically Symmetric in rest frame : general Spectator Particle (slow Proton)
      // **********************************************************
      //		psf     *= 4.*pi*pSMax^3/3;
      uw       = ran3.Uniform();
      pS_rest   = (pSMax)*pow(uw,1./3.); // uniform in 3p^2 dp = d(p^3)
      
      // Tel Aviv's high pT configuration
      //		pS_rest   = (pSMax+0.3)*pow(uu,1./3.); // uniform in 3p^2 dp = d(p^3)
      ux       = ran3.Uniform();
      csThRecoil = (2.*ux-1.);
      uy       = ran3.Uniform();
      phiRecoil  = pi*(2.*uy-1.);
      pS_3vec  = pS_rest*sin(acos(csThRecoil))*(cos(phiRecoil)*UnitXqCM + sin(phiRecoil)*UnitYqCM);
      pS_3vec += pS_rest*csThRecoil*UnitZqCM;
      pSpectator_Rest.SetVectM(pS_3vec,MSpectator);
      Jacob     *= 1./(2.*pSpectator_Rest.E());

      // pSpectator boost to Lab frame
      pSpectator_Vertex = pSpectator_Rest;
      pSpectator_Vertex.Boost(BoostRest);
		
      // "pSpectator_Rest" is the TLorentzVector of spectator particle in rest frame 
      // pSpectator_Rest.X(),pSpectator_Rest.Y(),pSpectator_Rest.Z(),pSpectator_Rest.E()
      //************************************************************************************
		
      PX_Vertex = kIncident_Vertex+PIncident_Vertex-kScattered_Vertex-pSpectator_Vertex;
      invts.MX2 = PX_Vertex.M2();

      // This spectator is not matter...
      invts.alphaS = ABeam*(pS_rest*csThRecoil+pSpectator_Rest.E())/MIon;
      invts.pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);
      invts.pPerpZ = pS_rest*sqrt(csThRecoil*csThRecoil);

      // For debugginh purpose
      /*
	cout << "(STEG variable #1) invariant variables: alphaS= " << invts.alphaS << endl;
      */
	    
      // alpha cut for select event with minimizing coherrent effects
      if(pow(invts.alphaS-1,2)<0.0001){
	continue;
      }

      // For debugging purpose
      /*
	cout  << "(6L) proton(TDIS spectator) Vertex, pRx= " << pprx_Lab << ", pPy= " << ppry_Lab <<
	", pPz= " << pprz_Lab << ", Epr= " << ppr_Lab << endl;
	cout  << "(6R) proton(TDIS spectator) Rest, pRx= " << pScatterProton_Rest.X() << ", pPy= " << pScatterProton_Rest.Y() <<
	", pPz= " << pScatterProton_Rest.Z() << ", Epr= " << pScatterProton_Rest.E() << endl;
	cout  << "(6L_S) proton(spectator) Vertex, pRx= " << pSpectator_Vertex.X()  << ", pPy= " << pSpectator_Vertex.Y()  <<
	", pPz= " << pSpectator_Vertex.Z()  << ", Epr= " << pSpectator_Vertex.E()  << endl;
	cout  << "(6R_S) proton(spectator) Rest, pRx= " << pSpectator_Rest.X()  << ", pPy= " << pSpectator_Rest.Y()  <<
	", pPz= " << pSpectator_Rest.Z()  << ", Epr= " << pSpectator_Rest.E()  << endl;
		
	cout << "Vertex: kScatteredPx=" << kScattered_Vertex.X() << ", kScatteredPy=" << kScattered_Vertex.Y() << ", kScatteredPz=" << kScattered_Vertex.Z() << endl;
	cout << "Rest: kScatteredPx=" << kScattered_Rest.X() << ", kScatteredPy=" << kScattered_Rest.Y() << ", kScatteredPz=" << kScattered_Rest.Z() << endl;
	      
	cout  << "(6L) pion(TDIS) Vertex, pRix= " << ppix_Lab << ", pPiy= " << ppiy_Lab <<
	", pPiz= " << ppiz_Lab << ", Epi= " << EpiE_Lab << endl;
	cout  << "(6R) pion(TDIS) Rest, pPix= " << Px_pi << ", pPiy= " << Py_pi <<
	", pPiz= " << Pz_pi << ", Epi= " << E_pi << endl;
      */
			      		
      invts.tSpectator = MIon*MIon+MSpectator*MSpectator-2.*pSpectator_Vertex.Dot(PIncident_Vertex);
      invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;
	    
      // This relation is from Ch. Weiss' note Eq.(35) , p_R**2
      invts.tempVar = (invts.tPrime/2)*(1+invts.tPrime/(2*MIon*MIon))+MIon*MIon/4-MSpectator*MSpectator;
	    
      //p_RT = sqrt((invts.tPrime)/2)*0.9; 
      invts.p_RT = sqrt((invts.tPrime)/2)*1.0; 

      //
      //  nominal Spectator Tagging Event Generator: either "proton" or "neutron" spectator particle
      //
      // store the 4-vector of the spectator nucleon in rest frame
      psx_Res = pSpectator_Rest.X() ;
      psy_Res = pSpectator_Rest.Y() ;
      psz_Res = pSpectator_Rest.Z() ;
      EsE_Res = pSpectator_Rest.E() ;
      vsx_Res = 0.0;
      vsy_Res = 0.0;
      vsz_Res = 0.0;
      
      // store the 4-vector of the spectator nucleon in Lab frame
      psx_Lab = pSpectator_Vertex.X() ;
      psy_Lab = pSpectator_Vertex.Y() ;
      psz_Lab = pSpectator_Vertex.Z() ;
      EsE_Lab = pSpectator_Vertex.E() ;
      vsx_Lab = 0.0;
      vsy_Lab = 0.0;
      vsz_Lab = 0.0;

    }
		
    PX_Vertex_Rest = kIncident_Rest+PIncident_Rest-kScattered_Rest-pSpectator_Rest;
	  
    // store the electron 4-vector in rest frame
    pex_Res = kScat3vec.X() ;
    pey_Res = kScat3vec.Y() ;
    pez_Res = kScat3vec.Z() ;
    EeE_Res = sqrt(kScat3vec.Mag2()+mElectron*mElectron);
    vex_Res = 0.0;
    vey_Res = 0.0;
    vez_Res = 0.0;
    
    // store the electron 4-vector in Lab frame
    pex_Lab = kScattered_Vertex.X() ;
    pey_Lab = kScattered_Vertex.Y() ;
    pez_Lab = kScattered_Vertex.Z() ;
    
    //EeE_Lab = sqrt(kScattered_Vertex.Mag2()+mElectron*mElectron);
    EeE_Lab = sqrt(pex_Lab*pex_Lab+pey_Lab*pey_Lab+pez_Lab*pez_Lab+mElectron*mElectron);
    vex_Lab = 0.0;
    vey_Lab = 0.0;
    vez_Lab = 0.0;

    // cout << "OUTPUT FOR GEMC " <<  Jacob << " \t"  <<  eBeamPol <<  " \t"  <<  DBeamPol << endl;	       

    if(ABeam>1.){
		  
      if( invts.alphaS<2.0 &&  invts.alphaS>0.0){
	OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  <<  NumPtls << " \t " <<  scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << "1.0" << " \t " << xpi << " \t" << ypi << " \t"  << tpi  << " \t"  <<  " \t" << sigma_dis << " \t" << sigma_tdis << endl;

	// incoming particles (e,D) = beam particles
	OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "1" << " \t " <<  e_particle_charge << " \t " << "0" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << kBeamMCx << " \t " << kBeamMCy << " \t " << kBeamMCz << " \t " << kBeamMC << " \t " << emass << " \t " << ve0X_Lab  << " \t " << ve0Y_Lab << " \t " << ve0Z_Lab << endl; 

	OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "2" << " \t " <<  d_particle_charge << " \t " << "0" << " \t " << d_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << PBeamMCx+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << PBeamMCy+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << PBeamMC*cos(CrossingTheta) << " \t " << PBeamMC << " \t " << MDeut << " \t " << vd0X_Lab  << " \t " << vd0Y_Lab << " \t " << vd0Z_Lab << endl; 

	// daughter particles
	// the scattered electron
	OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "3" << " \t " << e_particle_charge << " \t " << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific <<pex_Lab << " \t " << pey_Lab << " \t " << pez_Lab << " \t " << EeE_Lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
	// pion :  DIS: do NOT need to be detected
	OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "4" << " \t " << pi_particle_charge << " \t " << "1" << " \t " << pi_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific<< ppix_Lab << " \t " << ppiy_Lab << " \t " << ppiz_Lab << " \t " << EpiE_Lab << " \t " << mPion << " \t " << vpix_Lab  << " \t " << vpiy_Lab << " \t " << vpiz_Lab << endl; 
	// the first spectator protn (TDIS)
	OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "5" << " \t " << pr_particle_charge << " \t " << "1" << " \t " << pr_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific<< pprx_Lab << " \t " << ppry_Lab << " \t " << pprz_Lab << " \t " << EprE_Lab << " \t " << MProton << " \t " << vprx_Lab  << " \t " << vpry_Lab << " \t " << vprz_Lab << endl; 

      }
    }
    else{
      OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  <<  NumPtls << " \t " <<  scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << "1.0" << " \t " << xpi << " \t" << ypi << " \t"  << tpi  << " \t"  <<  " \t" << sigma_dis << " \t" << sigma_tdis << endl;
      
      // incoming particles (e,D) = beam particles
      OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "1" << " \t " << e_particle_charge << " \t " << "0" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << kBeamMCx << " \t " << kBeamMCy << " \t " << kBeamMCz << " \t " << kBeamMC << " \t " << emass << " \t " << ve0X_Lab  << " \t " << ve0Y_Lab << " \t " << ve0Z_Lab << endl; 

      OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "2" << " \t " << d_particle_charge << " \t " << "0" << " \t " << pr_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << PBeamMCx+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << PBeamMCy+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << PBeamMC*cos(CrossingTheta) << " \t " << PBeamMC << " \t " << MProton << " \t " << vd0X_Lab  << " \t " << vd0Y_Lab << " \t " << vd0Z_Lab << endl; 

      // daughter particles
      // the scattered electron
      OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "3" << " \t " << e_particle_charge << " \t " << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific <<pex_Lab << " \t " << pey_Lab << " \t " << pez_Lab << " \t " << EeE_Lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
      // pion :  DIS: do NOT need to be detected
      OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "4" << " \t " << pi_particle_charge << " \t " << "1" << " \t " << pi_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific<< ppix_Lab << " \t " << ppiy_Lab << " \t " << ppiz_Lab << " \t " << EpiE_Lab << " \t " << mPion << " \t " << vpix_Lab  << " \t " << vpiy_Lab << " \t " << vpiz_Lab << endl; 
      // the first spectator proton (TDIS)
      OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "5" << " \t " << pr_particle_charge << " \t " << "1" << " \t " << pr_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific<< pprx_Lab << " \t " << ppry_Lab << " \t " << pprz_Lab << " \t " << EprE_Lab << " \t " << MProton << " \t " << vprx_Lab  << " \t " << vpry_Lab << " \t " << vprz_Lab << endl; 

    }	
	  
    MEvts++;
    
    tree->Fill();

  }
  
  // itree->Fill();
  fRoot.Write();
  fRoot.Close();
  printf("Total of %d events out of %d Trials \n",MEvts,NEvts);
  
  OUT.close();

  return MEvts;
}
