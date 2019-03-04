int yields_kaon() {

  gROOT->Reset();

  TString DataDir = ".";

  //  TString myFile = DataDir + "/TTDIS-MC05x050.root"; // 5 on 50
  TString myFile1 = DataDir + "/TDIS-MC05x050_q23.75_kaon.root"; // 5 on 50, Q2=(3.5,4.0) GeV2
  TString myFile2 = DataDir + "/TDIS-MC05x050_q212_kaon.root"; // 5 on 50, Q2=(11.5,12.5) GeV2
  TString myFile3 = DataDir + "/TDIS-MC05x050_q230_kaon.root"; // 5 on 50, Q2=(28.8,31.2) GeV2
  TString myFile4 = DataDir + "/TDIS-MC05x050_q260_kaon.root"; // 5 on 50, Q2=(57.6,62.4) GeV2
  TString myFile5 = DataDir + "/TDIS-MC05x050_q2120_kaon.root"; // 5 on 50, Q2=(115,125) GeV2
  TString myFile6 = DataDir + "/TDIS-MC05x050_q2240_kaon.root"; // 5 on 50, Q2=(230,250) GeV2
 
  Char_t* PtrMyFile1 = (Char_t*)(myFile1);
  Char_t* PtrMyFile2 = (Char_t*)(myFile2);
  Char_t* PtrMyFile3 = (Char_t*)(myFile3);
  Char_t* PtrMyFile4 = (Char_t*)(myFile4);
  Char_t* PtrMyFile5 = (Char_t*)(myFile5);
  Char_t* PtrMyFile6 = (Char_t*)(myFile6);

  mySimFile1 = new TFile(PtrMyFile1);
  mySimFile2 = new TFile(PtrMyFile2);
  mySimFile3 = new TFile(PtrMyFile3);
  mySimFile4 = new TFile(PtrMyFile4);
  mySimFile5 = new TFile(PtrMyFile5);
  mySimFile6 = new TFile(PtrMyFile6);

  //Create an array for histograms (TObjArray)
  TObjArray *fHistoList = new TObjArray(50);

  c1 = new TCanvas("c1","GRAPH with ERRORS",200,10,700,500);
  //  TCanvas *c1 = new TCanvas("c1", "TDIS EIC MC",85,4,870,766);
  //TCanvas *c2 = new TCanvas("c2", "TDIS EIC MC",85,4,870,766);


  //create RECON histograms for real data                       

  TH1F* myXpiHist1 = new TH1F("myXpiHist1","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",10,0,1);
  TH1F* myXpiHist2 = new TH1F("myXpiHist2","5 on 50: xB=(0.001,1), Q2=12, s=1000",10,0,1);
  TH1F* myXpiHist3 = new TH1F("myXpiHist3","5 on 50: xB=(0.001,1), Q2=30, s=1000",10,0,1);
  TH1F* myXpiHist4 = new TH1F("myXpiHist4","5 on 50: xB=(0.001,1), Q2=60, s=1000",10,0,1);
  TH1F* myXpiHist5 = new TH1F("myXpiHist5","5 on 50: xB=(0.001,1), Q2=120, s=1000",10,0,1);
  TH1F* myXpiHist6 = new TH1F("myXpiHist6","5 on 50: xB=(0.001,1), Q2=240, s=1000",10,0,1);

  TH1F* myXpiHist7 = new TH1F("myXpiHist7","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",10,0,1);
  TH1F* myXpiHist8 = new TH1F("myXpiHist8","5 on 50: xB=(0.001,1), Q2=12, s=1000",10,0,1);
  TH1F* myXpiHist9 = new TH1F("myXpiHist9","5 on 50: xB=(0.001,1), Q2=30, s=1000",10,0,1);
  TH1F* myXpiHist10 = new TH1F("myXpiHist10","5 on 50: xB=(0.001,1), Q2=60, s=1000",10,0,1);
  TH1F* myXpiHist11 = new TH1F("myXpiHist11","5 on 50: xB=(0.001,1), Q2=120, s=1000",10,0,1);
  TH1F* myXpiHist12 = new TH1F("myXpiHist12","5 on 50: xB=(0.001,1), Q2=240, s=1000",10,0,1);

  TH1F* myfpi = new TH1F("myfpi","5 on 50: xB=(0.001,1), Q2=(1,100), s=1000",10,0,1);
  TH1F* myf2N = new TH1F("myf2N","5 on 50: xB=(0.001,1), Q2=(1,100), s=1000",10,0,1);


  //Cuts and Weights
  TCut cut1 = "sigma_dis";
  TCut cut2 = "sigma_tdis";

  TCut cut3 = "100000"; // total number of trials

  TCut totalCut = cut2 && cut3; // normalize MC output to cross section, where counts=integrated luminosity * cross setion. So the cross section represented by some event count is: count*total integrated cross section/total number of trials. The corresponding luminosity represented by a number of generated events is: total number of trials/total integrated cross section

  TCut cut11 = "xpi<0.1";
  TCut cut12 = "xpi>0.0";
  TCut xbinCut1 = cut11 && cut12;

  TCut cut13 = "xpi<0.2";
  TCut cut14 = "xpi>0.1";
  TCut xbinCut2 = cut13 && cut14;

  TCut cut15 = "xpi<0.3";
  TCut cut16 = "xpi>0.2";
  TCut xbinCut3 = cut15 && cut16;

  TCut cut17 = "xpi<0.4";
  TCut cut18 = "xpi>0.3";
  TCut xbinCut4 = cut17 && cut18;

  TCut cut19 = "xpi<0.5";
  TCut cut20 = "xpi>0.4";
  TCut xbinCut5 = cut19 && cut20;

  TCut cut21 = "xpi<0.6";
  TCut cut22 = "xpi>0.5";
  TCut xbinCut6 = cut21 && cut22;

  TCut cut23 = "xpi<0.7";
  TCut cut24 = "xpi>0.6";
  TCut xbinCut7 = cut23 && cut24;

  TCut cut25 = "xpi<0.8";
  TCut cut26 = "xpi>0.7";
  TCut xbinCut8 = cut25 && cut26;

  TCut cut27 = "xpi<0.9";
  TCut cut28 = "xpi>0.8";
  TCut xbinCut9 = cut27 && cut28;

  TCut cut29 = "xpi<1.0";
  TCut cut30 = "xpi>0.9";
  TCut xbinCut10 = cut29 && cut30;

  //split canvas                                                                          
  c1->Divide(3,2);

  //=======================================================================               

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c1->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Print();
  // myTree1->Draw("xpi>>myXpiHist1",totalCut);
  //myTree1->Draw("fpi>>myfpi");
  //myTree1->Draw("f2N>>myf2N");

  //Int_t counts=0;
  //counts =  myXpiHist1->GetEntries();
  //cout << "counts from 5 on 50: " << counts << endl;
  //fHistoList->AddAtFree(myXpiHist1);

  Int_t start=0;
  Int_t stop=0;
  Int_t incr=0;
  for (Int_t i=1; i<11; i=i+1){
    cout << "Index= " << i << endl;
    //    myTree1->Draw("xpi>>myXpiHist+i",totalCut&&xbinCut+i);
    //    myTree1->Draw("xpi>>myXpiHist+i",totalCut);
    myTree1->Draw("xpi>>myXpiHist1",totalCut);

  }

  myXpiHist1->GetXaxis()->SetTitle("xkaon");

  myXpiHist1->Draw();

 
  //myXpiHist->SetMarkerStyle(21);
  //myXpiHist->SetMarkerColor(kBlue);
  //myXpiHist->SetMarkerSize(0.5);
  //myXpiHist->SetLineColor(1);

  //myXpiHist->Draw("E1same");

  //c1->Update();


  //split canvas                                                                          
  //  c2->Divide(1,2);

  //c2->cd(1);
  //myfpi->Draw();
  //c2->cd(2);
  //myf2N->Draw();


  // XPI FOR 5 ON 50: Q2=12 GeV
  c1->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  myTree2->Draw("xpi>>myXpiHist2",totalCut);

  Int_t counts=0;
  counts =  myXpiHist2->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myXpiHist2);

  myXpiHist2->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c1->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  myTree3->Draw("xpi>>myXpiHist3",totalCut);

  Int_t counts=0;
  counts =  myXpiHist3->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myXpiHist3);

  myXpiHist3->Draw();



  // XPI FOR 5 ON 50: Q2=60 GeV
  c1->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  myTree4->Draw("xpi>>myXpiHist4",totalCut);

  Int_t counts=0;
  counts =  myXpiHist4->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myXpiHist4);

  myXpiHist4->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c1->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  myTree5->Draw("xpi>>myXpiHist5",totalCut);

  Int_t counts=0;
  counts =  myXpiHist5->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myXpiHist5);

  myXpiHist5->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c1->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  myTree6->Draw("xpi>>myXpiHist6",totalCut);

  Int_t counts=0;
  counts =  myXpiHist6->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myXpiHist6);

  myXpiHist6->Draw();

  c1->Print("kaon_distributions_5_on_50.pdf","pdf");

  gROOT->Reset();



  return 0;

}
