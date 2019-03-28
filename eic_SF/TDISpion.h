//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  1 16:49:17 2019 by ROOT version 6.12/06
// from TTree Evnts/p(e,e'np)X Event Generation   5 GeV/c x 100 GeV/c
// found on file: TDISpion.root
//////////////////////////////////////////////////////////

#ifndef TDISpion_h
#define TDISpion_h

#include <TROOT.h>
#include <TChain.h>
#include <TH2.h>
#include <TH3F.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TLorentzVector.h"



class TDISpion : public TSelector {
public :

  TString	foutname = "OUTPUT/plot_TDISpion";
  TString	outputpdf;

  double xbin[21]
    ={0.01,0.0125,0.0158,0.0199,0.0251,0.0316,0.0398,0.0501,0.0631,0.0794,0.1,0.125,0.158,0.199,0.251,0.316,0.398,0.501,0.631,0.794,1.00};
	
  double q2bin[11]
    ={1.0,1.58,2.51,3.98,6.31,10.0,15.8,25.1,39.8,63.1,100.};
  
  TTreeReader	 fReader;	//!the tree reader
  TTree	*fChain = 0;		//!pointer to the analyzed TTree or TChain
	
  TH1F	*inv_se;		//
  TH1F	*inv_sq;		//
  TH1F	*inv_Q2;		//  
  TH1F	*inv_xBj;		//
  TH1F	*inv_nu;		//
  TH1F	*inv_W;			//
  TH1F	*inv_x_D;		//
  TH1F	*inv_y_D;		//
  TH1F	*inv_tSpectator;	//
  TH1F	*inv_tPrime;		//		
  TH1F	*inv_TwoPdotk;		//
  TH1F	*inv_TwoPdotq;		//
  TH1F	*inv_RT;		//
  TH1F	*inv_pDrest;		//
  TH1F	*inv_tempVar;		//
  TH1F	*inv_MX2;		//
  TH1F	*inv_alphaS;		//
  TH1F	*inv_pPerpS;		//
  TH1F	*inv_pPerpZ;		//
	
  TH1F	*hxpi;			//
  TH1F	*hypi;			//
  TH1F	*htpi;			//
  TH1F	*hfpi;			//
  TH1F	*hf2N;			//
  TH1F	*hp2_pt;		//
  TH1F	*hp2_z;			//
  TH1F	*hsigma_dis;		//
  TH1F	*hsigma_tdis;		//
  TH1F	*hpex_Lab;		//
  TH1F	*hpey_Lab;		//
  TH1F	*hpez_Lab;		//
  TH1F	*hEeE_Lab;		//
  TH1F	*hppix_Lab;		//
  TH1F	*hppiy_Lab;		//
  TH1F	*hppiz_Lab;		//
  TH1F	*hEpiE_Lab;		//
  TH1F	*hpprx_Lab;		//
  TH1F	*hppry_Lab;		//
  TH1F	*hpprz_Lab;		//
  TH1F	*hEprE_Lab;		//
  TH1F	*hpXx_Lab;		//
  TH1F	*hpXy_Lab;		//
  TH1F	*hpXz_Lab;		//
  TH1F	*hEXE_Lab;		//
  TH1F	*hTDIS_xbj;		//
  TH1F	*hTDIS_znq;		//
  TH1F	*hTDIS_Mx2;		//
  TH1F	*hTDIS_y;		//
  TH1F	*hEScatRest;		//
  TH1F	*hkScatRest;		//
  TH1F	*hPhiScatRest;		//
  TH1F	*hcsPhiRest;		//
  TH1F	*hcsTheRest;		//

  TH1F *hbinned_tpi1;	//
  TH1F *hbinned_tpi2;	//
  TH1F *hbinned_tpi3;	//
  TH1F *hbinned_tpi4;	//

  TH2F *pex_v_csPhi;		//
  TH2F *pey_v_csPhi;		//
  TH2F *pez_v_csPhi;		//
  TH2F *q2_v_xbj;		//
  TH2F *f2N_v_q2;		//

  TH3F *hQ2vsXvsDISsigma; //
  TH3F *hQ2vsXvsDISsigma_cut; //
  

   // Readers to access the data (delete the ones you do not need).
   // TTreeReaderValue<TLorentzVector> e_Inc_ = {fReader, "e_Inc."};
   // TTreeReaderValue<TLorentzVector> P_Inc_ = {fReader, "P_Inc."};
   // TTreeReaderValue<TLorentzVector> e_Scat_ = {fReader, "e_Scat."};
   // TTreeReaderValue<TLorentzVector> p2_Sp_ = {fReader, "p2_Sp."};
   // TTreeReaderValue<TLorentzVector> q_Vir_ = {fReader, "q_Vir."};
   // TTreeReaderValue<TLorentzVector> p1_Sp_ = {fReader, "p1_Sp."};
   // TTreeReaderValue<TLorentzVector> pi_ = {fReader, "pi."};
   
   TTreeReaderValue<Double_t> s_e		= {fReader, "invts.s_e"};
   TTreeReaderValue<Double_t> s_q		= {fReader, "invts.s_q"};
   TTreeReaderValue<Double_t> Q2		= {fReader, "invts.Q2"};
   TTreeReaderValue<Double_t> xBj		= {fReader, "invts.xBj"};
   TTreeReaderValue<Double_t> nu		= {fReader, "invts.nu"};
   TTreeReaderValue<Double_t> W			= {fReader, "invts.W"};
   TTreeReaderValue<Double_t> x_D		= {fReader, "invts.x_D"};
   TTreeReaderValue<Double_t> y_D		= {fReader, "invts.y_D"};
   TTreeReaderValue<Double_t> tSpectator	= {fReader, "invts.tSpectator"};
   TTreeReaderValue<Double_t> tPrime		= {fReader, "invts.tPrime"};
   TTreeReaderValue<Double_t> TwoPdotk		= {fReader, "invts.TwoPdotk"};
   TTreeReaderValue<Double_t> TwoPdotq		= {fReader, "invts.TwoPdotq"};
   TTreeReaderValue<Double_t> p_RT		= {fReader, "invts.p_RT"};
   TTreeReaderValue<Double_t> pDrest		= {fReader, "invts.pDrest"};
   TTreeReaderValue<Double_t> tempVar		= {fReader, "invts.tempVar"};
   TTreeReaderValue<Double_t> MX2		= {fReader, "invts.MX2"};
   TTreeReaderValue<Double_t> alphaS		= {fReader, "invts.alphaS"};
   TTreeReaderValue<Double_t> pPerpS		= {fReader, "invts.pPerpS"};
   TTreeReaderValue<Double_t> pPerpZ		= {fReader, "invts.pPerpZ"};
   
   TTreeReaderValue<Double_t> xpi		= {fReader, "xpi"};
   TTreeReaderValue<Double_t> ypi		= {fReader, "ypi"};
   TTreeReaderValue<Double_t> tpi		= {fReader, "tpi"};
   TTreeReaderValue<Double_t> fpi		= {fReader, "fpi"};
   TTreeReaderValue<Double_t> f2N		= {fReader, "f2N"};
   TTreeReaderValue<Double_t> p2_pt		= {fReader, "p2_pt"};
   TTreeReaderValue<Double_t> p2_z		= {fReader, "p2_z"};
   TTreeReaderValue<Double_t> sigma_dis		= {fReader, "sigma_dis"};
   TTreeReaderValue<Double_t> sigma_tdis	= {fReader, "sigma_tdis"};
   TTreeReaderValue<Double_t> pex_Lab		= {fReader, "pex_Lab"};
   TTreeReaderValue<Double_t> pey_Lab		= {fReader, "pey_Lab"};
   TTreeReaderValue<Double_t> pez_Lab		= {fReader, "pez_Lab"};
   TTreeReaderValue<Double_t> EeE_Lab		= {fReader, "EeE_Lab"};
   TTreeReaderValue<Double_t> ppix_Lab		= {fReader, "ppix_Lab"};
   TTreeReaderValue<Double_t> ppiy_Lab		= {fReader, "ppiy_Lab"};
   TTreeReaderValue<Double_t> ppiz_Lab		= {fReader, "ppiz_Lab"};
   TTreeReaderValue<Double_t> EpiE_Lab		= {fReader, "EpiE_Lab"};
   TTreeReaderValue<Double_t> pprx_Lab		= {fReader, "pprx_Lab"};
   TTreeReaderValue<Double_t> ppry_Lab		= {fReader, "ppry_Lab"};
   TTreeReaderValue<Double_t> pprz_Lab		= {fReader, "pprz_Lab"};
   TTreeReaderValue<Double_t> EprE_Lab		= {fReader, "EprE_Lab"};
   TTreeReaderValue<Double_t> pXx_Lab		= {fReader, "pXx_Lab"};
   TTreeReaderValue<Double_t> pXy_Lab		= {fReader, "pXy_Lab"};
   TTreeReaderValue<Double_t> pXz_Lab		= {fReader, "pXz_Lab"};
   TTreeReaderValue<Double_t> EXE_Lab		= {fReader, "EXE_Lab"}; 
   TTreeReaderValue<Double_t> TDIS_xbj		= {fReader, "TDIS_xbj"};
   TTreeReaderValue<Double_t> TDIS_znq		= {fReader, "TDIS_znq"};
   TTreeReaderValue<Double_t> TDIS_Mx2		= {fReader, "TDIS_Mx2"};
   TTreeReaderValue<Double_t> TDIS_y		= {fReader, "TDIS_y"};
   TTreeReaderValue<Double_t> EScatRest		= {fReader, "EScatRest"};
   TTreeReaderValue<Double_t> kScatRest		= {fReader, "kScatRest"};
   TTreeReaderValue<Double_t> PhiScatRest	= {fReader, "PhiScatRest"};
   TTreeReaderValue<Double_t> csPhiRest		= {fReader, "csPhiRest"};
   TTreeReaderValue<Double_t> csTheRest		= {fReader, "csTheRest"};


   TDISpion(TTree * /*tree*/ =0) {inv_se=0,inv_sq=0,inv_Q2=0,  inv_xBj=0,inv_nu=0,inv_W=0,inv_x_D=0,inv_y_D=0,inv_tSpectator=0,inv_tPrime=0,inv_TwoPdotk=0,inv_TwoPdotq=0,inv_RT=0,inv_pDrest=0,inv_tempVar=0,inv_MX2=0,inv_alphaS=0,inv_pPerpS=0,inv_pPerpZ=0,hxpi=0,hypi=0,htpi=0,hfpi=0,hf2N=0,hp2_pt=0,hp2_z=0,hpex_Lab=0,hpey_Lab=0,hpez_Lab=0,hEeE_Lab=0,hppix_Lab=0,hppiy_Lab=0,hppiz_Lab=0,hEpiE_Lab=0,hpprx_Lab=0,hppry_Lab=0,hpprz_Lab=0,hEprE_Lab=0,hpXx_Lab=0,hpXy_Lab=0,hpXz_Lab=0,hEXE_Lab=0,hsigma_dis=0,hsigma_tdis=0,hTDIS_xbj=0,hTDIS_znq=0,hTDIS_Mx2=0,hTDIS_y=0,hEScatRest=0,hkScatRest=0,hPhiScatRest=0,hcsPhiRest=0,hcsTheRest=0,hbinned_tpi1=0,hbinned_tpi2=0,hbinned_tpi3=0,hbinned_tpi4=0,pex_v_csPhi=0,pey_v_csPhi=0,pez_v_csPhi=0,q2_v_xbj=0,f2N_v_q2=0,hQ2vsXvsDISsigma=0,hQ2vsXvsDISsigma_cut=0;}

   virtual ~TDISpion() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(TDISpion,0);

};

#endif

#ifdef TDISpion_cxx
void TDISpion::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t TDISpion::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef TDISpion_cxx
