#define TDISpion_cxx
// The class definition in TDISpion.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("TDISpion.C")
// root> T->Process("TDISpion.C","some options")
// root> T->Process("TDISpion.C+")
//


#include "TDISpion.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TText.h>

void TDISpion::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();
}

void TDISpion::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();

  inv_se		= new TH1F("inv_se","",200,1990,2010);						//
  inv_sq		= new TH1F("inv_sq","",200, 0,120);						//
  inv_Q2		= new TH1F("inv_Q2","",200,0,150);						//  
  inv_xBj		= new TH1F("inv_xBj","",200,0,1);						//
  inv_nu		= new TH1F("inv_nu","",200,0,100);						//
  inv_W			= new TH1F("inv_W","",200,-100,100);						//
  inv_x_D		= new TH1F("inv_x_D","",200,-1,1);						//
  inv_y_D		= new TH1F("inv_y_D","",200,-1,1);						//
  inv_tSpectator	= new TH1F("inv_tSpectator","",200,-1,1);					//
  inv_tPrime		= new TH1F("inv_tPrime","",200,0,3);						//
  inv_TwoPdotk		= new TH1F("inv_TwoPdotk","",200,0,0.1);					//
  inv_TwoPdotq		= new TH1F("inv_TwoPdotq","",200,-1,1);						//
  inv_RT		= new TH1F("inv_RT","",200,-1,1);						//
  inv_pDrest		= new TH1F("inv_pDrest","",200,1990,2010);					//
  inv_tempVar		= new TH1F("inv_tempVar","",200,0,120);						//
  inv_MX2		= new TH1F("inv_MX2","",200,-1,1);						//
  inv_alphaS		= new TH1F("inv_alphaS","",200,-1,1);						//
  inv_pPerpS		= new TH1F("inv_pPerpS","",200,-1,1);						//
  inv_pPerpZ		= new TH1F("inv_pPerpZ","",200,-1,1);						//
  
  hxpi			= new TH1F("hxpi","x_{#pi}",200,0,1.);						//
  hypi			= new TH1F("hypi","y_{#pi}",200,0,1.);						//
  htpi			= new TH1F("htpi","t_{#pi}",200,-1.,0);					        //
  hfpi			= new TH1F("hfpi","f_{#pi}",200,-1.,0);					        //
  hf2N			= new TH1F("hf2N","f2N",200,0,1.0);						//
  hp2_pt		= new TH1F("hp2_pt",".5% of incoming ion beam momentum (P_{iv})",200,-0.6,0.6);	//
  hp2_z			= new TH1F("hp2_z","Random number between (0,1)",200,-1.1,1.1);			//
  hsigma_dis		= new TH1F("hsigma_dis","DIS cross section",200,0,1000e3);			//
  hsigma_tdis		= new TH1F("hsigma_tdis","TDIS cross section",200,-200,200);			//
  hpex_Lab		= new TH1F("hpex_Lab","P_{ex} in Lab Frame",200,-15,15);			//
  hpey_Lab		= new TH1F("hpey_Lab","P_{ey} in Lab Frame",200,-15,15);			//
  hpez_Lab		= new TH1F("hpez_Lab","P_{ez} in Lab Frame",200,-15,15);			//
  hEeE_Lab		= new TH1F("hEeE_Lab","E_{e} in Lab Frame",200,0,15);				//
  hppix_Lab		= new TH1F("hppix_Lab","P_{#pix} in Lab Frame",200,-5,0);			//
  hppiy_Lab		= new TH1F("hppiy_Lab","P_{#piy} in Lab Frame",200,-0.6,0.6);			//
  hppiz_Lab		= new TH1F("hppiz_Lab","P_{#piz} in Lab Frame",200,-0.6,0.6);			//
  hEpiE_Lab		= new TH1F("hEpiE_Lab","E_{#pi} in Lab Frame",200,0,60);			//
  hpprx_Lab		= new TH1F("hpprx_Lab","P_{prx} in Lab Frame",200,-10,0);			//
  hppry_Lab		= new TH1F("hppry_Lab","P_{pry} in Lab Frame",200,-0.6,0.6);			//
  hpprz_Lab		= new TH1F("hpprz_Lab","P_{prz} in Lab Frame",200,-0.6,0.6);			//
  hEprE_Lab		= new TH1F("hEprE_Lab","E_{pr} in Lab Frame",200,0,200);			//
  hpXx_Lab		= new TH1F("hpXx_Lab","(Miss Mass) P_{Xx} in Lab Frame",200,-1,1);		//
  hpXy_Lab		= new TH1F("hpXy_Lab","(Miss Mass) P_{Xy} in Lab Frame",200,-1,1);		//
  hpXz_Lab		= new TH1F("hpXz_Lab","(Miss Mass) P_{Xz} in Lab Frame",200,-1,1);		//
  hEXE_Lab		= new TH1F("hEXE_Lab","(Miss Mass) E_{X} in Lab Frame",200,0,1);		//
  hTDIS_xbj		= new TH1F("hTDIS_xbj","TDIS x_{bj}",200,0,3);					//
  hTDIS_znq		= new TH1F("hTDIS_znq","",200,0,60);						//
  hTDIS_Mx2		= new TH1F("hTDIS_Mx2","",200,-100,100);					//
  hTDIS_y		= new TH1F("hTDIS_y","TDIS y",200,0,0.1);					//
  hEScatRest		= new TH1F("hEScatRest","Scat electron E_{e} in rest frame",200,900,1100);	//
  hkScatRest		= new TH1F("hkScatRest","Scat electron P_{e} in rest frame",200,900,1100);	//
  hPhiScatRest		= new TH1F("hPhiScatRest","Scat electron #phi in rest frame",200,-4,4);		//
  hcsPhiRest		= new TH1F("hcsPhiRest","Scat electron cos(#phi) in rest frame",200,-1.1,1.1);	//
  hcsTheRest		= new TH1F("hcsTheRest","Scat electron cos(#theta) in rest frame",200,0.9,1.1);	//

  hbinned_tpi1	= new TH1F("hbinned_tpi1","t binned in Q^{2} and x_{Bj}",200,0,3);		//
  hbinned_tpi2	= new TH1F("hbinned_tpi2","t binned in Q^{2} and x_{Bj}",200,0,3);		//
  hbinned_tpi3	= new TH1F("hbinned_tpi3","t binned in Q^{2} and x_{Bj}",200,0,3);		//
  hbinned_tpi4	= new TH1F("hbinned_tpi4","t binned in Q^{2} and x_{Bj}",200,0,3);		//

  pex_v_csPhi		= new TH2F("pex_v_csPhi","P_{ex} vs. #phi",200,-4,4,200,-15,15);		//
  pey_v_csPhi		= new TH2F("pey_v_csPhi","P_{ey} vs. #phi",200,-4,4,200,-15,15);		//
  pez_v_csPhi		= new TH2F("pez_v_csPhi","P_{ez} vs. #phi",200,-4,4,200,-15,15);		//
  q2_v_xbj		= new TH2F("q2_v_xbj","Q^{2} vs. x_{Bj}",200,0,1,200,0,150);		//
  f2N_v_q2		= new TH2F("f2pi_v_q2","f2N vs. Q^{2}",200,0,150,200,0,1.0);		//
  

  GetOutputList()->Add(inv_se);			//
  GetOutputList()->Add(inv_sq);			//
  GetOutputList()->Add(inv_Q2);			//  
  GetOutputList()->Add(inv_xBj);		//
  GetOutputList()->Add(inv_nu);			//
  GetOutputList()->Add(inv_W);			//
  GetOutputList()->Add(inv_x_D);		//
  GetOutputList()->Add(inv_y_D);		//
  GetOutputList()->Add(inv_tSpectator);		//
  GetOutputList()->Add(inv_tPrime);		//  		
  GetOutputList()->Add(inv_TwoPdotk);		//
  GetOutputList()->Add(inv_TwoPdotq);		//
  GetOutputList()->Add(inv_RT);			//
  GetOutputList()->Add(inv_pDrest);		//
  GetOutputList()->Add(inv_tempVar);		//
  GetOutputList()->Add(inv_MX2);		//
  GetOutputList()->Add(inv_alphaS);		//
  GetOutputList()->Add(inv_pPerpS);		//
  GetOutputList()->Add(inv_pPerpZ);		//
  
  GetOutputList()->Add(hxpi);			//
  GetOutputList()->Add(hypi);			//
  GetOutputList()->Add(htpi);			//
  GetOutputList()->Add(hfpi);			//
  GetOutputList()->Add(hf2N);			//
  GetOutputList()->Add(hp2_pt);			//
  GetOutputList()->Add(hp2_z);			//
  GetOutputList()->Add(hsigma_dis);		//
  GetOutputList()->Add(hsigma_tdis);		//
  GetOutputList()->Add(hpex_Lab);		//
  GetOutputList()->Add(hpey_Lab);		//
  GetOutputList()->Add(hpez_Lab);		//
  GetOutputList()->Add(hEeE_Lab);		//
  GetOutputList()->Add(hppix_Lab);		//
  GetOutputList()->Add(hppiy_Lab);		//
  GetOutputList()->Add(hppiz_Lab);		//
  GetOutputList()->Add(hEpiE_Lab);		//
  GetOutputList()->Add(hpprx_Lab);		//
  GetOutputList()->Add(hppry_Lab);		//
  GetOutputList()->Add(hpprz_Lab);		//
  GetOutputList()->Add(hEprE_Lab);		//
  GetOutputList()->Add(hpXx_Lab);		//
  GetOutputList()->Add(hpXy_Lab);		//
  GetOutputList()->Add(hpXz_Lab);		//
  GetOutputList()->Add(hEXE_Lab);		//
  GetOutputList()->Add(hTDIS_xbj);		//
  GetOutputList()->Add(hTDIS_znq);		//
  GetOutputList()->Add(hTDIS_Mx2);		//
  GetOutputList()->Add(hTDIS_y);		//
  GetOutputList()->Add(hEScatRest);		//
  GetOutputList()->Add(hkScatRest);		//
  GetOutputList()->Add(hPhiScatRest);		//
  GetOutputList()->Add(hcsPhiRest);		//
  GetOutputList()->Add(hcsTheRest);		//
  
  GetOutputList()->Add(hbinned_tpi1);	//
  GetOutputList()->Add(hbinned_tpi2);	//
  GetOutputList()->Add(hbinned_tpi3);	//
  GetOutputList()->Add(hbinned_tpi4);	//

  GetOutputList()->Add(pex_v_csPhi);		//
  GetOutputList()->Add(pey_v_csPhi);		//
  GetOutputList()->Add(pez_v_csPhi);		//
  GetOutputList()->Add(q2_v_xbj);		//
  GetOutputList()->Add(f2N_v_q2);		//
}

Bool_t TDISpion::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
   
   inv_se->Fill(*s_e);				//
   inv_sq->Fill(*s_q);				//
   inv_Q2->Fill(*Q2);				//  
   inv_xBj->Fill(*xBj);				//
   inv_nu->Fill(*nu);				//
   inv_W->Fill(*W);				//
   inv_x_D->Fill(*x_D);				//
   inv_y_D->Fill(*y_D);				//
   inv_tSpectator->Fill(*tSpectator);		//
   inv_tPrime->Fill(*tPrime);			//
   inv_TwoPdotk->Fill(*TwoPdotk);		//
   inv_TwoPdotq->Fill(*TwoPdotq);		//
   inv_RT->Fill(*p_RT);				//
   inv_pDrest->Fill(*pDrest);			//
   inv_tempVar->Fill(*tempVar);			//
   inv_MX2->Fill(*MX2);				//
   inv_alphaS->Fill(*alphaS);			//
   inv_pPerpS->Fill(*pPerpS);			//
   inv_pPerpZ->Fill(*pPerpZ);			//
  
   hxpi->Fill(*xpi);				//
   hypi->Fill(*ypi);				//
   htpi->Fill(*tpi);				//
   hfpi->Fill(*fpi);				//
   hf2N->Fill(*f2N);				//
   hp2_pt->Fill(*p2_pt);			//
   hp2_z->Fill(*p2_z);				//
   hsigma_dis->Fill(*sigma_dis);		//
   hsigma_tdis->Fill(*sigma_tdis);		//
   hpex_Lab->Fill(*pex_Lab);			//
   hpey_Lab->Fill(*pey_Lab);			//
   hpez_Lab->Fill(*pez_Lab);			//
   hEeE_Lab->Fill(*EeE_Lab);			//
   hppix_Lab->Fill(*ppix_Lab);			//
   hppiy_Lab->Fill(*ppiy_Lab);			//
   hppiz_Lab->Fill(*ppiz_Lab);			//
   hEpiE_Lab->Fill(*EpiE_Lab);			//
   hpprx_Lab->Fill(*pprx_Lab);			//
   hppry_Lab->Fill(*ppry_Lab);			//
   hpprz_Lab->Fill(*pprz_Lab);			//
   hEprE_Lab->Fill(*EprE_Lab);			//
   hpXx_Lab->Fill(*pXx_Lab);			//
   hpXy_Lab->Fill(*pXy_Lab);			//
   hpXz_Lab->Fill(*pXz_Lab);			//
   hEXE_Lab->Fill(*EXE_Lab);			//
   hTDIS_xbj->Fill(*TDIS_xbj);			//
   hTDIS_znq->Fill(*TDIS_znq);			//
   hTDIS_Mx2->Fill(*TDIS_Mx2);			//
   hTDIS_y->Fill(*TDIS_y);			//
   hEScatRest->Fill(*EScatRest);		//
   hkScatRest->Fill(*kScatRest);		//
   hPhiScatRest->Fill(*PhiScatRest);		//
   hcsPhiRest->Fill(*csPhiRest);		//
   hcsTheRest->Fill(*csTheRest);		//

   pex_v_csPhi->Fill(*PhiScatRest,*pex_Lab);	//
   pey_v_csPhi->Fill(*PhiScatRest,*pey_Lab);	//
   pez_v_csPhi->Fill(*PhiScatRest,*pez_Lab);	//
   q2_v_xbj->Fill(*xBj,*Q2,1.);	//
   f2N_v_q2->Fill(*Q2,*f2N);   //
   
   if((*xBj >= 0. && *xBj <= 1.0) &&
      (*Q2 >= 0 && *Q2 <= 5)
      ) {
     
     hbinned_tpi4->Fill(*tpi);		//
    
   }

   if((*xBj >= 0. && *xBj <= 1.0) &&
      (*Q2 >= 60 && *Q2 <= 65)
      ) {
     
     hbinned_tpi3->Fill(*tpi);		//
    
   }

   if((*xBj >= 0. && *xBj <= 0.5) &&
      (*Q2 >= 0 && *Q2 <= 100)
      ) {
     
     hbinned_tpi2->Fill(*tpi);		//
    
   }

   if((*xBj >= 0.5 && *xBj <= 1.) &&
      (*Q2 >= 0 && *Q2 <= 100)
      ) {
     
     hbinned_tpi1->Fill(*tpi);		//
    
   }
   
   return kTRUE;
}

void TDISpion::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TDISpion::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  outputpdf = foutname + ".pdf";

  TCanvas		*c_inv_se;	        //
  c_inv_se		= new TCanvas("c_inv_se","");
  c_inv_se->cd();
  inv_se->Draw();
  c_inv_se->Print(outputpdf + "(");
  TCanvas		*c_inv_sq;		//
  c_inv_sq		= new TCanvas("c_inv_sq","");
  c_inv_sq->cd();
  inv_sq->Draw();
  c_inv_sq->Print(outputpdf);
  TCanvas		*c_inv_Q2;		//
  c_inv_Q2		= new TCanvas("c_inv_Q2","");
  c_inv_Q2->cd();
  inv_Q2->Draw();
  c_inv_Q2->Print(outputpdf);
  TCanvas		*c_inv_xBj;		//
  c_inv_xBj		= new TCanvas("c_inv_xBj","");
  c_inv_xBj->cd();
  inv_xBj->Draw();
  c_inv_xBj->Print(outputpdf);
  TCanvas		*c_inv_nu;		//
  c_inv_nu		= new TCanvas("c_inv_nu","");
  c_inv_nu->cd();
  inv_nu->Draw();
  c_inv_nu->Print(outputpdf);
  TCanvas		*c_inv_W;		//
  c_inv_W		= new TCanvas("c_inv_W","");
  c_inv_W->cd();
  inv_W->Draw();
  c_inv_W->Print(outputpdf);
  TCanvas		*c_inv_x_D;		//
  c_inv_x_D		= new TCanvas("c_inv_x_D","");
  c_inv_x_D->SetLogy();
  c_inv_x_D->cd();
  inv_x_D->Draw();
  inv_x_D->Print(outputpdf);
  TCanvas		*c_inv_y_D;		//
  c_inv_y_D		= new TCanvas("c_inv_y_D","");
  c_inv_y_D->SetLogy();
  c_inv_y_D->cd();
  inv_y_D->Draw();
  inv_y_D->Print(outputpdf);
  TCanvas		*c_inv_tSpectator;	//
  c_inv_tSpectator	= new TCanvas("c_inv_tSpectator","");
  c_inv_tSpectator->SetLogy();
  c_inv_tSpectator->cd();
  inv_tSpectator->Draw();
  c_inv_tSpectator->Print(outputpdf);
  TCanvas		*c_inv_tPrime;		//
  c_inv_tPrime		= new TCanvas("c_inv_tPrime","");
  c_inv_tPrime->cd();
  inv_tPrime->Draw();
  c_inv_tPrime->Print(outputpdf);
  TCanvas		*c_inv_TwoPdotk;	//
  c_inv_TwoPdotk	= new TCanvas("c_inv_TwoPdotk","");
  c_inv_TwoPdotk->cd();
  inv_TwoPdotk->Draw();
  c_inv_TwoPdotk->Print(outputpdf);
  TCanvas		*c_inv_TwoPdotq;	//
  c_inv_TwoPdotq	= new TCanvas("c_inv_TwoPdotq","");
  c_inv_TwoPdotq->SetLogy();
  c_inv_TwoPdotq->cd();
  inv_TwoPdotq->Draw();
  c_inv_TwoPdotq->Print(outputpdf);
  TCanvas		*c_inv_RT;		//
  c_inv_RT		= new TCanvas("c_inv_RT","");
  c_inv_RT->SetLogy();
  c_inv_RT->cd();
  inv_RT->Draw();
  c_inv_RT->Print(outputpdf);
  TCanvas		*c_inv_pDrest;		//
  c_inv_pDrest		= new TCanvas("c_inv_pDrest","");
  c_inv_pDrest->cd();
  inv_pDrest->Draw();
  c_inv_pDrest->Print(outputpdf);
  TCanvas		*c_inv_tempVar;		//
  c_inv_tempVar		= new TCanvas("c_inv_tempVar","");
  c_inv_tempVar->cd();
  inv_tempVar->Draw();
  c_inv_tempVar->Print(outputpdf);
  TCanvas		*c_inv_MX2;		//
  c_inv_MX2		= new TCanvas("c_inv_MX2","");
  c_inv_MX2->SetLogy();
  c_inv_MX2->cd();
  inv_MX2->Draw();
  c_inv_MX2->Print(outputpdf);
  TCanvas		*c_inv_alphaS;		//
  c_inv_alphaS		= new TCanvas("c_inv_alphaS","");
  c_inv_alphaS->SetLogy();
  c_inv_alphaS->cd();
  inv_alphaS->Draw();
  c_inv_alphaS->Print(outputpdf);
  TCanvas		*c_inv_pPerpS;		//
  c_inv_pPerpS		= new TCanvas("c_inv_pPerpS","");
  c_inv_pPerpS->SetLogy();
  c_inv_pPerpS->cd();
  inv_pPerpS->Draw();
  c_inv_pPerpS->Print(outputpdf);
  TCanvas		*c_inv_pPerpZ;		//
  c_inv_pPerpZ		= new TCanvas("c_inv_pPerpZ","");
  c_inv_pPerpZ->SetLogy();
  c_inv_pPerpZ->cd();
  inv_pPerpZ->Draw();
  c_inv_pPerpZ->Print(outputpdf);
  
  TCanvas		*c_hxpi;		//
  c_hxpi		= new TCanvas("c_hxpi","");
  c_hxpi->SetLogy();
  c_hxpi->cd();
  hxpi->Draw();
  c_hxpi->Print(outputpdf);
  TCanvas		*c_hypi;		//
  c_hypi		= new TCanvas("c_hypi","");
  c_hypi->cd();
  hypi->Draw();
  c_hypi->Print(outputpdf);
  TCanvas		*c_htpi;		//
  c_htpi		= new TCanvas("c_htpi","");
  c_htpi->SetLogy();
  c_htpi->cd();
  htpi->Draw();
  c_htpi->Print(outputpdf);
  TCanvas		*c_hfpi;		//
  c_hfpi		= new TCanvas("c_hfpi","");
  c_hfpi->SetLogy();
  c_hfpi->cd();
  hfpi->Draw();
  c_hfpi->Print(outputpdf);
  TCanvas		*c_hf2N;		//
  c_hf2N		= new TCanvas("c_hf2N","");
  c_hf2N->cd();
  hf2N->Draw();
  c_hf2N->Print(outputpdf);
  TCanvas		*c_hp2_pt;		//
  c_hp2_pt		= new TCanvas("c_hp2_pt","");
  c_hp2_pt->cd();
  hp2_pt->Draw();
  c_hp2_pt->Print(outputpdf);
  TCanvas		*c_hp2_z;		//
  c_hp2_z		= new TCanvas("c_hp2_z","");
  c_hp2_z->cd();
  hp2_z->Draw();
  c_hp2_z->Print(outputpdf);
  TCanvas		*c_hsigma_dis;		//
  c_hsigma_dis		= new TCanvas("c_hsigma_dis","");
  c_hsigma_dis->SetLogy();
  c_hsigma_dis->cd();
  hsigma_dis->Draw();
  c_hsigma_dis->Print(outputpdf);
  TCanvas		*c_hsigma_tdis;		//
  c_hsigma_tdis		= new TCanvas("c_hsigma_tdis","");
  c_hsigma_tdis->SetLogy();
  c_hsigma_tdis->cd();
  hsigma_tdis->Draw();
  c_hsigma_tdis->Print(outputpdf);
  TCanvas		*c_hpex_Lab;		//
  c_hpex_Lab		= new TCanvas("c_hpex_Lab","");
  //c_hpex_Lab->SetLogy();
  c_hpex_Lab->cd();
  hpex_Lab->Draw();
  c_hpex_Lab->Print(outputpdf);
  TCanvas		*c_hpey_Lab;		//
  c_hpey_Lab		= new TCanvas("c_hpey_Lab","");
  //c_hpey_Lab->SetLogy();
  c_hpey_Lab->cd();
  hpey_Lab->Draw();
  c_hpey_Lab->Print(outputpdf);
  TCanvas		*c_hpez_Lab;		//
  c_hpez_Lab		= new TCanvas("c_hpez_Lab","");
  //c_hpez_Lab->SetLogy();
  c_hpez_Lab->cd();
  hpez_Lab->Draw();
  c_hpez_Lab->Print(outputpdf);
  TCanvas		*c_hEeE_Lab;		//
  c_hEeE_Lab		= new TCanvas("c_hEeE_Lab","");
  //c_hEeE_Lab->SetLogy();
  c_hEeE_Lab->cd();
  hEeE_Lab->Draw();
  c_hEeE_Lab->Print(outputpdf);
  TCanvas		*c_hppix_Lab;		//
  c_hppix_Lab		= new TCanvas("c_hppix_Lab","");
  //c_hppix_Lab->SetLogy();
  c_hppix_Lab->cd();
  hppix_Lab->Draw();
  c_hppix_Lab->Print(outputpdf);
  TCanvas		*c_hppiy_Lab;		//
  c_hppiy_Lab		= new TCanvas("c_hppiy_Lab","");
  //c_hppiy_Lab->SetLogy();
  c_hppiy_Lab->cd();
  hppiy_Lab->Draw();
  c_hppiy_Lab->Print(outputpdf);
  TCanvas		*c_hppiz_Lab;		//
  c_hppiz_Lab		= new TCanvas("c_hppiz_Lab","");
  //c_hppiz_Lab->SetLogy();
  c_hppiz_Lab->cd();
  hppiz_Lab->Draw();
  c_hppiz_Lab->Print(outputpdf);
  TCanvas		*c_hEpiE_Lab;		//
  c_hEpiE_Lab		= new TCanvas("c_hEpiE_Lab","");
  //c_hEpiE_Lab->SetLogy();
  c_hEpiE_Lab->cd();
  hEpiE_Lab->Draw();
  c_hEpiE_Lab->Print(outputpdf);
  TCanvas		*c_hpprx_Lab;		//
  c_hpprx_Lab		= new TCanvas("c_hpprx_Lab","");
  c_hpprx_Lab->SetLogy();
  c_hpprx_Lab->cd();
  hpprx_Lab->Draw();
  c_hpprx_Lab->Print(outputpdf);
  TCanvas		*c_hppry_Lab;		//
  c_hppry_Lab		= new TCanvas("c_hppry_Lab","");
  c_hppry_Lab->SetLogy();
  c_hppry_Lab->cd();
  hppry_Lab->Draw();
  c_hppry_Lab->Print(outputpdf);
  TCanvas		*c_hpprz_Lab;		//
  c_hpprz_Lab		= new TCanvas("c_hpprz_Lab","");
  c_hpprz_Lab->SetLogy();
  c_hpprz_Lab->cd();
  hpprz_Lab->Draw();
  c_hpprz_Lab->Print(outputpdf);
  TCanvas		*c_hEprE_Lab;		//
  c_hEprE_Lab		= new TCanvas("c_hEprE_Lab","");
  c_hEprE_Lab->SetLogy();
  c_hEprE_Lab->cd();
  hEprE_Lab->Draw();
  c_hEprE_Lab->Print(outputpdf);
  TCanvas		*c_hpXx_Lab;		//
  c_hpXx_Lab		= new TCanvas("c_hpXx_Lab","");
  c_hpXx_Lab->SetLogy();
  c_hpXx_Lab->cd();
  hpXx_Lab->Draw();
  c_hpXx_Lab->Print(outputpdf);
  TCanvas		*c_hpXy_Lab;		//
  c_hpXy_Lab		= new TCanvas("c_hpXy_Lab","");
  c_hpXy_Lab->SetLogy();
  c_hpXy_Lab->cd();
  hpXy_Lab->Draw();
  c_hpXy_Lab->Print(outputpdf);
  TCanvas		*c_hpXz_Lab;		//
  c_hpXz_Lab		= new TCanvas("c_hpXz_Lab","");
  c_hpXz_Lab->SetLogy();
  c_hpXz_Lab->cd();
  hpXz_Lab->Draw();
  c_hpXz_Lab->Print(outputpdf);
  TCanvas		*c_hEXE_Lab;		//
  c_hEXE_Lab		= new TCanvas("c_hEXE_Lab","");
  c_hEXE_Lab->SetLogy();
  c_hEXE_Lab->cd();
  hEXE_Lab->Draw();
  c_hEXE_Lab->Print(outputpdf);
  TCanvas		*c_hTDIS_xbj;		//
  c_hTDIS_xbj		= new TCanvas("c_hTDIS_xbj","");
  c_hTDIS_xbj->cd();
  hTDIS_xbj->Draw();
  c_hTDIS_xbj->Print(outputpdf);
  TCanvas		*c_hTDIS_znq;		//
  c_hTDIS_znq		= new TCanvas("c_hTDIS_znq","");
  c_hTDIS_znq->cd();
  hTDIS_znq->Draw();
  c_hTDIS_znq->Print(outputpdf);
  TCanvas		*c_hTDIS_Mx2;		//
  c_hTDIS_Mx2		= new TCanvas("c_hTDIS_Mx2","");
  c_hTDIS_Mx2->cd();
  hTDIS_Mx2->Draw();
  c_hTDIS_Mx2->Print(outputpdf);
  TCanvas		*c_hTDIS_y;		//
  c_hTDIS_y		= new TCanvas("c_hTDIS_y","");
  c_hTDIS_y->cd();
  hTDIS_y->Draw();
  c_hTDIS_y->Print(outputpdf);
  TCanvas		*c_EScatRest;		//
  c_EScatRest		= new TCanvas("c_EScatRest","");
  c_EScatRest->cd();
  hEScatRest->Draw();
  c_EScatRest->Print(outputpdf);
  TCanvas		*c_kScatRest;		//
  c_kScatRest		= new TCanvas("c_kScatRest","");
  c_kScatRest->cd();
  hkScatRest->Draw();
  c_kScatRest->Print(outputpdf);
  TCanvas		*c_PhiScatRest;		//
  c_PhiScatRest		= new TCanvas("c_PhiScatRest","");
  c_PhiScatRest->cd();
  hPhiScatRest->Draw();
  c_PhiScatRest->Print(outputpdf);
  TCanvas		*c_csPhiRest;		//
  c_csPhiRest		= new TCanvas("c_csPhiRest","");
  c_csPhiRest->cd();
  hcsPhiRest->Draw();
  c_csPhiRest->Print(outputpdf);
  TCanvas		*c_csTheRest;		//
  c_csTheRest		= new TCanvas("c_csTheRest","");
  c_csTheRest->cd();
  hcsTheRest->Draw();
  c_csTheRest->Print(outputpdf);

  
  ////////////////////////////////////////////////
  TCanvas		*c_hbinned_tpi;    	//
  TPaveText *tab; //
  c_hbinned_tpi		= new TCanvas("c_hbinned_tpi","t binned in Q^{2} and x_{Bj}");
  tab = new TPaveText(.95,.45,.3,.65,"NDC");
  c_hbinned_tpi->SetLogy();
  c_hbinned_tpi->cd();
  gStyle->SetOptStat(0);
  hbinned_tpi1->SetLineColor(2);
  hbinned_tpi1->SetFillColor(46);
  hbinned_tpi2->SetLineColor(3);
  hbinned_tpi2->SetFillColor(30);
  hbinned_tpi3->SetLineColor(4);
  hbinned_tpi3->SetFillColor(38);
  hbinned_tpi4->SetLineColor(6);
  hbinned_tpi4->SetFillColor(47);
  hbinned_tpi2->Draw();
  hbinned_tpi1->Draw("SAME");
  hbinned_tpi3->Draw("SAME");
  hbinned_tpi4->Draw("SAME");
  c_hbinned_tpi->Update();
  TText *l1 = tab->AddText("(*xBj >= 0.5 && *xBj <= 1.) && (*Q2 >= 0 && *Q2 <= 100)"); l1->SetTextColor(2);
  TText *l2 = tab->AddText("(*xBj >= 0. && *xBj <= 0.5) && (*Q2 >= 0 && *Q2 <= 100)"); l2->SetTextColor(3);
  TText *l3 = tab->AddText("(*xBj >= 0. && *xBj <= 1.0) && (*Q2 >= 60 && *Q2 <= 65)"); l3->SetTextColor(4);
  TText *l4 = tab->AddText("(*xBj >= 0. && *xBj <= 1.0) && (*Q2 >= 0 && *Q2 <= 5)"); l4->SetTextColor(6);
  tab->Draw();
  c_hbinned_tpi->Print(outputpdf);
									    
  TCanvas		*c_pex_v_csPhi;		//
  c_pex_v_csPhi		= new TCanvas("c_pex_v_csPhi","");
  c_pex_v_csPhi->cd();
  pex_v_csPhi->Draw();
  c_pex_v_csPhi->Print(outputpdf);

  TCanvas		*c_pey_v_csPhi;		//
  c_pey_v_csPhi		= new TCanvas("c_pey_v_csPhi","");
  c_pey_v_csPhi->cd();
  pey_v_csPhi->Draw();
  c_pey_v_csPhi->Print(outputpdf);

  TCanvas		*c_pez_v_csPhi;		//
  c_pez_v_csPhi		= new TCanvas("c_pez_v_csPhi","");
  c_pez_v_csPhi->cd();
  pez_v_csPhi->Draw();
  c_pez_v_csPhi->Print(outputpdf);

  TCanvas		*c_q2_v_xbj;		//
  c_q2_v_xbj		= new TCanvas("c_q2_v_xbj","");
  c_q2_v_xbj->cd();
  q2_v_xbj->Draw("colz");
  c_q2_v_xbj->Print(outputpdf);

  TCanvas		*c_f2N_v_q2;		//
  c_f2N_v_q2		= new TCanvas("c_f2N_v_q2","");
  c_f2N_v_q2->cd();
  f2N_v_q2->Draw("colz");
  c_f2N_v_q2->Print(outputpdf+")");

}
