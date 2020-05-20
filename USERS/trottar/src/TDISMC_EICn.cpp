/*
 * Description: Electron on Proton beam, tagged neutron and pi+ final states 
 *              Please see README for instructions
 * ================================================================
 * Time-stamp: "2020-05-20 15:22:21 trottar"
 * ================================================================
 *
 * Author:  Kijun Park and Richard L. Trotta III <trotta@cua.edu>
 *
 * Copyright (c) trottar
 */

// General definition in this header
#include "TDISMC_EIC.h"
#include "structure_functions/tim_hobbs/timhobbs.h"
#include "structure_functions/sf.h"
#include <time.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <TROOT.h>

// init DIS cteq pdf
void initcteqpdf();

int getIndexOfLargestElement(Double_t arr[], int FromIndex, int ToIndex);

// Define form factor function for proton and pion : f2p and f2pi to get TDIS cross-section
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

  // Define the DIS PDF from CTEQ directory:  cteq-tbls/ctq66m/ctq66.00.pds
  initcteqpdf();
  
  TRandom3 ran3;// Random number generator [45ns/call]
  ran3.SetSeed(rnum);// Sets random generator seed which is used initialize the random number generator
  
  Bool_t iran=kTRUE;      // TRUE==include incident beam emittance
  //Bool_t iran=kFALSE;      // Falut NO incident beam emittance

  // Define particle properties
  double emass =  mElectron, prmass = MProton, MSpectator = MNeutron,
    spmass = MSpectator, mPion = mPip, pimass = mPion;
  int e_particle_charge = -1, pr_particle_charge = 1, sp_particle_charge = 0, pi_particle_charge = 1;
  int e_particle_id = 11, pr_particle_id =  2212, sp_particle_id = 2112, pi_particle_id = 211;

  // Number of incident and final state particles of reaction
  int NumPtls = 5;

  // 1 for proton, 2 for neutron incident
  int inucl = 1;

  // Incident proton information
  double pprx_inc,ppry_inc,pprz_inc,EprE_inc ;

  // vertex for the initial electron
  double_t ve0X_Lab,ve0Y_Lab,ve0Z_Lab ;
  
  // vertex for the initial proton
  double_t vp0X_Lab,vp0Y_Lab,vp0Z_Lab ;

  // TDIS scattered electron information
  double pex_Res,pey_Res,pez_Res,EeE_Res;
  double pex_Lab,pey_Lab,pez_Lab,EeE_Lab ;
  double_t vex_Lab,vey_Lab,vez_Lab ;
  double_t vex_Res,vey_Res,vez_Res ;
 
  // TDIS detected pion information
  double ppix_Res,ppiy_Res,ppiz_Res,EpiE_Res;
  double ppix_Lab,ppiy_Lab,ppiz_Lab,EpiE_Lab ;
  double_t vpix_Lab,vpiy_Lab,vpiz_Lab ;
  double_t vpix_Res,vpiy_Res,vpiz_Res ;
 
  // TDIS spectator neutron information
  double pnx_Lab,pny_Lab,pnz_Lab,EnE_Lab ;
  double_t vnx_Lab,vny_Lab,vnz_Lab ;

  // missing particle (X) information
  double pXx_Res,pXy_Res,pXz_Res,EXE_Res ;
  double_t vXx_Res,vXy_Res,vXz_Res ;
  double pXx_Lab,pXy_Lab,pXz_Lab,EXE_Lab ;
  double_t vXx_Lab,vXy_Lab,vXz_Lab ;

  TLorentzVector kIncident_Vertex, PIncident_Vertex, qVirtual_Vertex, kScattered_Vertex;

  // Define daughter particle Vertex
  TLorentzVector pSpectator_Vertex,  PX_Vertex, PX_Vertex_Rest;
  TLorentzVector pScatterNeutron_Vertex;
  TLorentzVector pScatterPion_Vertex;

  // TDIS-SBS xbj means xbj_D (xbj of nucleon should be xbj/2 ? Asking Dasuni for xbj definition
  double Px_p2,Py_p2,Pz_p2,P_p2,theta_p2,phi_p2,E_p2, p2_pt, p2_z;

  double xVpiR, xSpiR;
  double yy,rho_PS, W2;
  double F2pi,kTk,pi_int,F2piK;
  
  // Pion structure function selection
  //  typ = 1 ! s-dependent expon FORM FACTOR (L=1.56,dis=1,FLAG=0)
  //  typ = 2 ! COV DIP FORM FACTOR
  //  typ = 3 ! t-dependent expon FORM FACTOR (L=0.85,dis=1,FLAG=0)
  //  typ = 4 ! Paul-Villas FORM FACTOR (L=0.27,dis=1,FLAG=0)
  //  typ = 5 ! ZEUS parameterization with F2N (proton SF)
  double typ = 5;

  double weight_tdis;
  
  char tTitle[80], tName[18], rName[32];
  
  sprintf(tTitle,"p(e,e'\u03C0n)X Event Generation %3.0f GeV/c x%4.0f GeV/c",kBeam, PBeam);
  sprintf(rName,"../OUTPUTS/pi_n_%.0fon%.0f.root", kBeam, PBeam);

  TFile fRoot(rName,"Recreate", tTitle);
  sprintf(tName,"pi_n");

  TTree *tree = new TTree("Evnts",tTitle);

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
  tree->Branch("p1_Sp.", &pScatterNeutron_Vertex, bufsize, 1);
  tree->Branch("pi.", &pScatterPion_Vertex, bufsize, 1);

  // Store physics varaibles into ROOTs
  tree->Branch("xpi", &xpi, "xpi/D");
  tree->Branch("ypi", &ypi, "ypi/D");
  tree->Branch("y", &y, "y/D");
  tree->Branch("t", &y, "t/D");
  tree->Branch("tpi", &tpi, "tpi/D");
  tree->Branch("fpi", &fpi, "fpi/D");
	
  tree->Branch("f2N", &f2N, "f2N/D");
	
  tree->Branch("sigma_dis", &sigma_dis, "sigma_dis/D");
  tree->Branch("sigma_tdis", &sigma_tdis, "sigma_tdis/D");

  // Incident proton
  tree->Branch("pprx_inc", &pprx_inc, "pprx_inc/D");
  tree->Branch("ppry_inc", &ppry_inc, "ppry_inc/D");
  tree->Branch("pprz_inc", &pprz_inc, "pprz_inc/D");
  tree->Branch("EprE_inc", &EprE_inc, "EprE_inc/D");
  
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

  // TDIS spectator neutron information
  tree->Branch("pnx_Lab", &pnx_Lab, "pnx_Lab/D");
  tree->Branch("pny_Lab", &pny_Lab, "pny_Lab/D");
  tree->Branch("pnz_Lab", &pnz_Lab, "pnz_Lab/D");
  tree->Branch("EnE_Lab", &EnE_Lab, "EnE_Lab/D");

  // missing particle (X) information
  tree->Branch("pXx_Lab", &pXx_Lab, "pXx_Lab/D");
  tree->Branch("pXy_Lab", &pXy_Lab, "pXy_Lab/D");
  tree->Branch("pXz_Lab", &pXz_Lab, "pXz_Lab/D");
  tree->Branch("EXE_Lab", &EXE_Lab, "EXE_Lab/D");

  // TDIS kinematic variables
  tree->Branch("TDIS_xbj", &TDIS_xbj, "TDIS_xbj/D");
  tree->Branch("TDIS_znq", &TDIS_znq, "TDIS_znq/D");
  tree->Branch("TDIS_Mx2", &TDIS_Mx2, "TDIS_Mx2/D");
  tree->Branch("TDIS_y", &TDIS_y, "TDIS_y/D");
	
  double sig_eThx=0.0, sig_eThy=0.0, sig_iThx=0.0, sig_iThy=0.0;

  double MIon = MProton;
  
  // Global Phase Space Factor;
  //s_e - M_ion^2 - m_electron^2
  double uu, uv, uw, ux, uy, norm;
  
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
  PIncident_0.SetXYZM(PIon*sin(CrossingTheta)*cos(CrossingPhi),
		      PIon*sin(CrossingTheta)*sin(CrossingPhi),PIon*cos(CrossingTheta),MIon);
  
  TTree *itree = new TTree("Meta",tTitle);
  itree->Branch("e0.", &kIncident_0,bufsize,1);
  itree->Branch("P0.",  &PIncident_0, bufsize, 1);
  
  typedef struct {
    Int_t nEvts;
    Float_t PhSpFct;
  } MonteCarlo;
  static	MonteCarlo mc;
  mc.nEvts = NEvts;
  mc.PhSpFct   = (Q2Max-Q2Min)*(log(xMax)-log(xMin))*4.*pi*pow((pSMax),3)/3.;
  itree->Branch("MC", &mc, "nEvts/I:PhSpFct/F");

  typedef struct {
    Double_t p2_pt, p2_z,
      phi_p2, theta_p2,
      Pz_p2, P_p2, E_p2,
      Px_p2, Py_p2,Pz_p2_1,
      Pz_p2_2,Pz_p2_3,Pz_p2_4; 
  } P2;
  static P2 p2;

  itree->Branch("P2", &p2, "p2_pt/D:p2_z/D:phi_p2/D:theta_p2/D:Pz_p2/D:P_p2/D:E_p2/D:Px_p2/D:Py_p2/D:Pz_p2_1/D:Pz_p2_2/D:Pz_p2_3/D:Pz_p2_4/D");
  itree->Branch("Jacob",&Jacob,"Jacob/D");

  double pS_rest, csThRecoil, phiRecoil;

  //name of output lund file for GEANT4 use
  ofstream OUT (Form("../OUTPUTS/pi_n_%.0fon%.0f_lund.dat", kBeam, PBeam), ios::trunc);

  // Get LorentzVector for the pScattered Proton for TDIS in rest frame
  TLorentzVector pScatterNeutron_Rest;

  double pn_Lab;

  // *****************************************************************************
  // Generate missing Hadron(Pion) for TDIS in rest frame
  // initial neutron 4momentum(-P1) - TDIS scattered proton 4momentum(P2)
  // *****************************************************************************
  double E_pi,Px_pi,Py_pi,Pz_pi;

  TLorentzVector pScatterPion_Rest;

  TVector3 pScatterPion_V3, pScatterNeutron_V3;

  double P_pi;

  double fpifac = 1;

  double minE = 0.0;
  
  //for progress bar
  double progress=0.0;

  Int_t MEvts=1;
  
  // Start of event generator
  for (int iEvt=0; iEvt<=NEvts; iEvt++) {

    // Progress bar
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
    vp0X_Lab = 0;
    vp0Y_Lab = 0;
    vp0Z_Lab = 0;	 	
		
    // these are with beam smearing
    kIncident_Vertex.SetXYZM(kBeamMCx, kBeamMCy, kBeamMCz, mElectron);// Set 4-momentum of incident electron beam
    PIncident_Vertex.SetXYZM(PBeamMCx, PBeamMCy, PBeamMCz, MIon);// Set 4-momentum of incident ion beam
    
    //  Crossing angle ehric with smearing
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
    
    // EprE_inc = PIncident_Vertex.E();
    EprE_inc = sqrt(pprx_inc*pprx_inc+ppry_inc*ppry_inc+pprz_inc*pprz_inc+MProton*MProton);

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
    uu   = ran3.Uniform();
    invts.Q2   = Q2Max*uu + Q2Min*(1.-uu);
    uv   = ran3.Uniform();
    invts.xBj  = pow(xMin,1.-uv)*pow(xMax,uv);
    invts.x_D  = invts.xBj*(MProton/MIon);
    invts.y_D  = invts.Q2/(invts.x_D*invts.TwoPdotk);
    invts.Yplus = 1 + ((1-invts.y_D)*(1-invts.y_D));

    y = invts.y_D;

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    if (invts.y_D>=(1.0-2.*mElectron*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
    if (invts.y_D>=1./(1.+invts.x_D*MIon*MIon/invts.TwoPdotk) ) {
      // Unphysical kinematics
      continue;
    }
	  
    Jacob *=invts.xBj;   // Jacobian  dx/dlnx

    // scattered electron, phi angle in rest frame 
    PhiScatRest = pi*(2.*ran3.Uniform()-1.0);
    csPhiRest = cos(PhiScatRest);
    //  mElectron-->0 approximation
    EScatRest = invts.TwoPdotk*(1.-invts.y_D)/(2.*MIon);

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
    if (EScatRest<mElectron) {
      // should never happen
      printf("illegal Rest frame scattered electron energy =%10.6f \n",EScatRest);
      continue;
    }
    
    kScatRest = sqrt(EScatRest*EScatRest-mElectron*mElectron);// scattered electron momentum in rest frame
    csTheRest = (2.*EScatRest*kIncident_Rest.E() - invts.Q2 - 2.*mElectron*mElectron) /(2.*kScatRest*kIncident_Rest.P());   // scattered electron \cosine\theta in rest frame	                                                 	  

    // **** for debugging purpose
    // for the dubugging purpose: OK, CHECKED
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
	  
    // Proton momentum in rest frame with and without smearing
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

    // Null for proton beam
    pSpectator_Rest  = PIncident_Rest;
    
    //  definition are moved at the beginning of code
    //		double TDIS_xbj, TDIS_znq,TDIS_Mx2,TDIS_y;
    TDIS_xbj = invts.Q2/(2*pSpectator_Rest.Dot(qVirtual_Rest));
    TDIS_znq = p2.p2_z*pSpectator_Rest.Dot(qVirtual_Rest);
    TDIS_Mx2 = (qVirtual_Rest + pSpectator_Rest).Mag2();
    TDIS_y   = (pSpectator_Rest.Dot(qVirtual_Rest))/(pSpectator_Rest.Dot(kIncident_Rest));

    //***
    // For debugging purpose:
    /*
      cout << "(3) kinematics: TDIS.xBj = " << TDIS_xbj <<  ", invts.xBj= " << invts.xBj <<  ", invts.Q2= " << invts.Q2<< endl; 
    */
    
    // Not needed!
    // neutron config = proton + pion(-):  this proton called additional speactator in TDIS concept
    // randomize in pt and z
    p2_pt = gRandom->Uniform(0.005*PBeam); // .5% of incoming ion beam momentum, this is the limit of the transverse momentum of  recoil particle...
    p2_z = gRandom->Uniform(1.);
    phi_p2 = gRandom->Uniform(360.0*D2R);
    Px_p2 = p2_pt*cos(phi_p2);
    Py_p2 = p2_pt*sin(phi_p2);
    Pz_p2 = ( -1.0*TDIS_znq*qVirtual_Rest.Z() + sqrt( TDIS_znq*qVirtual_Rest.Z()*TDIS_znq*qVirtual_Rest.Z() + invts.Q2*(qVirtual_Rest.E()*qVirtual_Rest.E())*(MProton*MProton + p2_pt*p2_pt)- invts.Q2*TDIS_znq*TDIS_znq) )/invts.Q2/1000.;  // unit match with "p2_pt" (GeV)
    P_p2  = sqrt (Pz_p2*Pz_p2 + p2_pt*p2_pt);
    theta_p2 = acos (Pz_p2/P_p2);
    E_p2 = sqrt ( P_p2*P_p2 + MProton*MProton);
  
    // Define scattered proton kinematics
    // Generate pSpectator Spherically Symmetric in rest frame     
    uw       = ran3.Uniform();
    pS_rest   = (pSMax)*pow(uw,1./3.); // uniform in 3p^2 dp = d(p^3), pSMax=0.3
    ux       = ran3.Uniform();
    csThRecoil = (2.*ux-1.);
    uy       = ran3.Uniform();
    phiRecoil  = pi*(2.*uy-1.);
    
    pScatterNeutron_V3  = pS_rest*sin(acos(csThRecoil))*(cos(phiRecoil)*UnitXqCM
							+ sin(phiRecoil)*UnitYqCM);
    pScatterNeutron_V3 += pS_rest*csThRecoil*UnitZqCM;
    
    pScatterNeutron_Rest.SetVectM(pScatterNeutron_V3,MSpectator);
    Jacob     *= 1./(2.*pScatterNeutron_Rest.E());

    PX_Vertex_Rest = kIncident_Rest+PIncident_Rest-(kScattered_Rest+pScatterNeutron_Rest);

    invts.MX2 = PX_Vertex.M2();
    invts.alphaS = ABeam*(pS_rest*csThRecoil+pSpectator_Rest.E())/MIon;
    invts.pPerpS = pS_rest*sqrt(1.-csThRecoil*csThRecoil);

    invts.tSpectator = MIon*MIon+MSpectator*MSpectator - 2.*pSpectator_Vertex.Dot(PIncident_Vertex);
    invts.tPrime     = 2.*pSpectator_Vertex.Dot(PIncident_Vertex) - MIon*MIon;

    // alpha cut for select event with minimizing coherrent effects
    if(pow(invts.alphaS-1,2)<0.0001){
      continue;
    }

	  
    // for the debugging purpose: SEEMS NOT CRAZY NUMBER....
    // if (TMath::Sqrt(TDIS_Mx2) > PBeam/2){
    //   cout << "---->TDIS missing mass =" << TMath::Sqrt(TDIS_Mx2)  << endl;
    // }
      
    E_pi  = pSpectator_Rest.E() - pScatterNeutron_Rest.E();
    Px_pi = pSpectator_Rest.X() - pScatterNeutron_Rest.X(); 
    Py_pi = pSpectator_Rest.Y() - pScatterNeutron_Rest.Y();
    Pz_pi = pSpectator_Rest.Z() - pScatterNeutron_Rest.Z();
      
    // pScatterPion_Rest.SetXYZM(Px_pi,Py_pi,Pz_pi,mPion);
    pScatterPion_Rest = pSpectator_Rest - (pScatterNeutron_Rest+PX_Vertex_Rest);
    // pScatterPion_Rest = pSpectator_Rest - pScatterNeutron_Rest;
      
    pScatterPion_V3.SetXYZ(Px_pi,Py_pi,Pz_pi);

    // Back to Lab frame
    pScatterNeutron_Vertex = pScatterNeutron_Rest;	  
    pScatterNeutron_Vertex.Boost(BoostRest);   
    pScatterPion_Vertex = pScatterPion_Rest;
    pScatterPion_Vertex.Boost(BoostRest);

    P_pi = pScatterPion_V3.Mag();

    if ((pScatterNeutron_Vertex.E()+pScatterPion_Vertex.E()) > PBeam){
    // if ((pScatterNeutron_Vertex.E()) > PBeam){
      continue;
    }
    
    pnx_Lab = pScatterNeutron_Vertex.X();
    pny_Lab = pScatterNeutron_Vertex.Y();
    pnz_Lab = pScatterNeutron_Vertex.Z();
    EnE_Lab = sqrt(pnx_Lab*pnx_Lab+pny_Lab*pny_Lab+pnz_Lab*pnz_Lab+MProton*MProton);
      
    vnx_Lab = 0.0;
    vny_Lab = 0.0;
    vnz_Lab = 0.0;
    
    pn_Lab = sqrt(pnx_Lab*pnx_Lab+pny_Lab*pny_Lab+pnz_Lab*pnz_Lab);
      
    ppix_Lab = pScatterPion_Vertex.X();
    ppiy_Lab = pScatterPion_Vertex.Y();
    ppiz_Lab = pScatterPion_Vertex.Z();
    EpiE_Lab = sqrt(ppix_Lab*ppix_Lab+ppiy_Lab*ppiy_Lab+ppiz_Lab*ppiz_Lab+mPion*mPion);
    
    vpix_Lab = 0.0;
    vpiy_Lab = 0.0;
    vpiz_Lab = 0.0;
        
    // For debugging purpose
    // if (pn_Lab > (PBeam/ABeam) ) { 
    //   // Unphysical kinematics  // if TDIS spectator momentum larger than 50% ion momentum
    //   printf("impossible of TDIS spectator momentum= %6.2f \n", ppr_Lab);
    //   continue;
    // })

    // for debuggin purpose
    // cout  << "(7L) pion Vertex, px= " << ppix_Lab << ", py= " << ppiy_Lab <<
    // ", pz= " << ppiz_Lab << ", Epi= " << EpiE_Lab << endl;
	       
    // Definition of pion variable
    xpi = TDIS_xbj/(1 - p2.p2_z);
    tpi = (E_pi*E_pi) - (pScatterPion_V3.Mag()*pScatterPion_V3.Mag());
    ypi = pScatterPion_Rest.Dot(qVirtual_Rest)/(pScatterPion_Rest.Dot(kIncident_Rest));
    
    // *****
    // for debugging purpose
    /*
      cout  << "(7R) pion REST, px= " << pScatterPion_Rest.X() << ", py= " << pScatterPion_Rest.Y()<<
      ", pz= " << pScatterPion_Rest.Z() << ", E= " << pScatterPion_Rest.E() << ", M= " << pScatterPion_Rest.M() << endl;
    */
	  
    // if (TDIS_xbj > 0.001 && TDIS_xbj < 1.0){
    if (xpi > 0.001 && xpi < 1.0){ // xpi for ZEUS parameterization
      
      /*
	L = renormalization cut-off parameter... ex) typ = 1 ! 1.56GeV
	dis = charge state configuration (0=charge exchange, 1=neutral exchange)
	FLAG = combination of splitting function/PDF (0:\pi, 1:\rho, 2:\pi\Delta)
      */
      //  typ = 1 ! s-dependent expon FORM FACTOR (L=1.56,dis=1,FLAG=0)      
      if (typ == 1){
      fpi = fpifac*f2pi(P_pi, TDIS_xbj, theta_p2/D2R);
      }
      
      //  typ = 2 ! COV DIP FORM FACTOR      
      if (typ == 2){
      fpi = fpifac*f2picov(P_pi, TDIS_xbj, theta_p2/D2R);
      }
      
      //  typ = 3 ! t-dependent expon FORM FACTOR (L=0.85,dis=1,FLAG=0)      
      if (typ == 3){
      fpi = fpifac*f2pitexp(P_pi, TDIS_xbj,theta_p2/D2R);
      }
      
      //  typ = 4 ! Paul-Villas FORM FACTOR (L=0.27,dis=1,FLAG=0)      
      if (typ == 4){
      fpi = fpifac*f2pipaulvillas(P_pi, TDIS_xbj, theta_p2/D2R);
      }
      
      //  typ = 5 ! ZEUS parameterization with F2N (proton SF)      
      if (typ == 5){
	fpi = fpifac*f2piZEUS(TDIS_xbj,invts.Q2,inucl);
      }

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

      f2N   = F2N(invts.xBj, invts.Q2, inucl);  // F2N is define at G4SBSDIS.h
      sigma_dis    = cdissigma_p(TDIS_xbj,invts.y_D,invts.Q2,invts.nu,kScattered_Rest.E());
      // sigma_dis    = cdissigma_p(TDIS_xbj,TDIS_y,invts.Q2,invts.nu,kScattered_Rest.E());
      sigma_tdis   = sigma_dis * (fpi/f2N); 
    }

    // ***
    // For debugging purpose
    /*
      if(fpi>0.){
      cout << " what is f2N= "  << f2N  << ", fpi=" << fpi << ", momentum = " << P_pi << ", theta = " << theta_p2/D2R << endl;
      }
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
		
    PX_Vertex = kIncident_Vertex+PIncident_Vertex-kScattered_Vertex-pSpectator_Vertex;
	  
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
    
    EeE_Lab = sqrt(pex_Lab*pex_Lab+pey_Lab*pey_Lab+pez_Lab*pez_Lab+mElectron*mElectron);
    vex_Lab = 0.0;
    vey_Lab = 0.0;
    vez_Lab = 0.0;

    OUT << setiosflags(ios::left)  << setiosflags(ios::fixed)  <<"                 "  <<  NumPtls << " \t " <<  scientific  << invts.xBj << " \t " << invts.Q2  << " \t " << invts.s_e  << " \t " << "1.0" << " \t " << xpi << " \t" << ypi << " \t"  << tpi  << " \t"  <<  " \t" << sigma_dis << " \t" << sigma_tdis << endl;
      
    // incoming particles
    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "1" << " \t " << e_particle_charge << " \t " << "0" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << kBeamMCx << " \t " << kBeamMCy << " \t " << kBeamMCz << " \t " << kBeamMC << " \t " << emass << " \t " << ve0X_Lab  << " \t " << ve0Y_Lab << " \t " << ve0Z_Lab << endl; 

    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "2" << " \t " << pr_particle_charge << " \t " << "0" << " \t " << pr_particle_id << " \t " << "0" <<  " \t "<< "1" << " \t "<< scientific << PBeamMCx+PBeamMC*sin(CrossingTheta)*cos(CrossingPhi) << " \t " << PBeamMCy+PBeamMC*sin(CrossingTheta)*sin(CrossingPhi) << " \t " << PBeamMC*cos(CrossingTheta) << " \t " << PBeamMC << " \t " << MProton << " \t " << vp0X_Lab  << " \t " << vp0Y_Lab << " \t " << vp0Z_Lab << endl; 

    // daughter particles
    // the scattered electron
    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) << "\t" << "3" << " \t " << e_particle_charge << " \t " << "1" << " \t " << e_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << pex_Lab << " \t " << pey_Lab << " \t " << pez_Lab << " \t " << EeE_Lab << " \t " << emass << " \t " << vex_Lab  << " \t " << vey_Lab << " \t " << vez_Lab << endl; 
    // final state pion
    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "4" << " \t " << pi_particle_charge << " \t " << "1" << " \t " << pi_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << ppix_Lab << " \t " << ppiy_Lab << " \t " << ppiz_Lab << " \t " << EpiE_Lab << " \t " << mPion << " \t " << vpix_Lab  << " \t " << vpiy_Lab << " \t " << vpiz_Lab << endl; 
    // final state neutron (TDIS)
    OUT << setiosflags(ios::left) << setiosflags(ios::fixed) <<  "\t" << "5" << " \t " << sp_particle_charge << " \t " << "1" << " \t " << sp_particle_id << " \t " << "0" <<  " \t "<< "1" <<  " \t "<< scientific << pnx_Lab << " \t " << pny_Lab << " \t " << pnz_Lab << " \t " << EnE_Lab << " \t " << spmass << " \t " << vnx_Lab  << " \t " << vny_Lab << " \t " << vnz_Lab << endl;  
    
    MEvts++;
    
    tree->Fill();
    itree->Fill();

  }
  
  fRoot.Write();
  fRoot.Close();
  printf("Total of %d events out of %d Trials \n",MEvts,NEvts);
  
  OUT.close();

  return MEvts;
}
