int yields_kaon_xbin() {

  gROOT->Reset();

  TString DataDir = ".";

  TString myFile11 = DataDir + "/TDIS-MC05x050_q23.75_kaon.root"; // 5 on 50, Q2=(3.6,3.9) GeV2
  TString myFile21 = DataDir + "/TDIS-MC05x050_q212_kaon.root"; // 5 on 50, Q2=(11.5,12.5) GeV2
  TString myFile31 = DataDir + "/TDIS-MC05x050_q230_kaon.root"; // 5 on 50, Q2=(28.7,31.3) GeV2
  TString myFile41 = DataDir + "/TDIS-MC05x050_q260_kaon.root"; // 5 on 50, Q2=(57.5,62.5) GeV2
  TString myFile51 = DataDir + "/TDIS-MC05x050_q2120_kaon.root"; // 5 on 50, Q2=(115,125) GeV2
  TString myFile61 = DataDir + "/TDIS-MC05x050_q2240_kaon.root"; // 5 on 50, Q2=(230,250) GeV2

  TString myFile1 = DataDir + "/TDIS-MC05x100_q23.75_kaon.root"; // 5 on 100, Q2=(3.6,3.9) GeV2
  TString myFile2 = DataDir + "/TDIS-MC05x100_q212_kaon.root"; // 5 on 100, Q2=(11.5,12.5) GeV2
  TString myFile3 = DataDir + "/TDIS-MC05x100_q230_kaon.root"; // 5 on 100, Q2=(28.7,31.3) GeV2
  TString myFile4 = DataDir + "/TDIS-MC05x100_q260_kaon.root"; // 5 on 100, Q2=(57.5,62.5) GeV2
  TString myFile5 = DataDir + "/TDIS-MC05x100_q2120_kaon.root"; // 5 on 100, Q2=(115,125) GeV2
  TString myFile6 = DataDir + "/TDIS-MC05x100_q2240_kaon.root"; // 5 on 100, Q2=(230,250) GeV2
 
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
  c2 = new TCanvas("c2","GRAPH with ERRORS 2",0,100,700,500);
  c3 = new TCanvas("c3","GRAPH with ERRORS 3",100,200,700,500);
  c4 = new TCanvas("c4","GRAPH with ERRORS 4",200,300,700,500);
  c5 = new TCanvas("c5","GRAPH with ERRORS 5",300,400,700,500);
  c6 = new TCanvas("c6","GRAPH with ERRORS 6",400,500,700,500);
  c7 = new TCanvas("c7","GRAPH with ERRORS 7",0,100,800,600);
  c8 = new TCanvas("c8","GRAPH with ERRORS 8",0,100,900,700);
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

  TH1F* myfpi1 = new TH1F("myfpi1","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",100,-0.0001,0.0001);
  TH1F* myfpi2 = new TH1F("myfpi2","5 on 50: xB=(0.001,1), Q2=12, s=1000",100,-0.0001,0.0001);
  TH1F* myfpi3 = new TH1F("myfpi3","5 on 50: xB=(0.001,1), Q2=30, s=1000",100,-0.0001,0.0001);
  TH1F* myfpi4 = new TH1F("myfpi4","5 on 50: xB=(0.001,1), Q2=60, s=1000",100,-0.0001,0.0001);
  TH1F* myfpi5 = new TH1F("myfpi5","5 on 50: xB=(0.001,1), Q2=120, s=1000",100,-0.0001,0.0001);
  TH1F* myfpi6 = new TH1F("myfpi6","5 on 50: xB=(0.001,1), Q2=240, s=1000",100,-0.0001,0.0001);

  TH1F* myfpix21 = new TH1F("myfpix21","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix22 = new TH1F("myfpix22","5 on 50: xB=(0.001,1), Q2=12, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix23 = new TH1F("myfpix23","5 on 50: xB=(0.001,1), Q2=30, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix24 = new TH1F("myfpix24","5 on 50: xB=(0.001,1), Q2=60, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix25 = new TH1F("myfpix25","5 on 50: xB=(0.001,1), Q2=120, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix26 = new TH1F("myfpix26","5 on 50: xB=(0.001,1), Q2=240, s=1000",100,-0.0001,0.0001);

  TH1F* myfpix31 = new TH1F("myfpix31","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix32 = new TH1F("myfpix32","5 on 50: xB=(0.001,1), Q2=12, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix33 = new TH1F("myfpix33","5 on 50: xB=(0.001,1), Q2=30, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix34 = new TH1F("myfpix34","5 on 50: xB=(0.001,1), Q2=60, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix35 = new TH1F("myfpix35","5 on 50: xB=(0.001,1), Q2=120, s=1000",100,-0.0001,0.0001);
  TH1F* myfpix36 = new TH1F("myfpix36","5 on 50: xB=(0.001,1), Q2=240, s=1000",100,-0.0001,0.0001);

  TH1F* myf2N1 = new TH1F("myf2N1","5 on 50: xB=(0.001,1), Q2=3.75, s=1000",10,-1,2);


  //Cuts and Weights
  TCut cut1 = "sigma_dis";
  TCut cut2 = "sigma_tdis";

  TCut cut3 = "2.86E-05"; // inverse of total number of trials (=35000)
  TCut cut4 = "8.6E+10"; // Luminosity
  TCut cut5 = "1E-06"; // factor to normalize to number of events in bin rather than overall number (cut3)

  TCut totalCut = cut2 && cut3 && cut4; // normalize MC output to cross section, where counts=integrated luminosity * cross section. So the cross section represented by some event count is: count*total integrated cross section/total number of trials. The corresponding luminosity represented by a number of generated events is: total number of trials/total integrated cross section


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
  c2->Divide(3,2);
  c3->Divide(3,2);
  c4->Divide(3,2);
  c5->Divide(3,2);
  c6->Divide(3,2);
  c7->Divide(3,2);
  c8->Divide(3,2);

  //=======================================================================               

  //XBIN 1

  cout << "======= XBIN 1: x=0.05 =====" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c1->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
    myTree1->Print();
    //myTree1->Draw("fpi>>myfpi1",totalCut&&xbinCut1);
    myTree1->Draw("fpi>>myfpi1",xbinCut1);

    cout << "fpi mean for Q2=3.75, x=0.05 is: " << myfpi1.GetMean() << endl;
    //cout << "fpi RMS for Q2=3.75 is: " << myfpi1.GetRMS() << endl;

  Int_t counts=0;
  counts =  myfpi1->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi1);

  myXpiHist1->GetXaxis()->SetTitle("xkaon");

 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c1->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpi2",totalCut&&xbinCut1);
  myTree2->Draw("fpi>>myfpi2",xbinCut1);

    cout << "fpi mean for Q2=12, x=0.05 is: " << myfpi2.GetMean() << endl;
    //cout << "fpi RMS for Q2=12 is: " << myfpi2.GetRMS() << endl;

  Int_t counts=0;
  counts =  myfpi2->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi2);

  // myfpi2->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c1->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpi3",totalCut&&xbinCut1);  
   myTree3->Draw("fpi>>myfpi3",xbinCut1); 

    cout << "fpi mean for Q2=30, x=0.05 is: " << myfpi3.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpi3->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi3);

  //myfpi3->Draw();



  // XPI FOR 5 ON 50: Q2=60 GeV
  c1->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  //myTree4->Draw("fpi>>myfpi4",totalCut&&xbinCut1); 
  myTree4->Draw("fpi>>myfpi4",xbinCut1);

    cout << "fpi mean for Q2=60, x=0.05 is: " << myfpi4.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpi4->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi4);

  //myfpi4->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c1->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  // myTree5->Draw("fpi>>myfpi5",totalCut&&xbinCut1);
    myTree5->Draw("fpi>>myfpi5",xbinCut1);

    cout << "fpi mean for Q2=120, x=0.05 is: " << myfpi5.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpi5->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi5);

  //myfpi5->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c1->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpi6",totalCut&&xbinCut1);
    myTree6->Draw("fpi>>myfpi6",xbinCut1);

  cout << "fpi mean for Q2=240, x=0.05 is: " << myfpi6.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpi6->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpi6);

  //myfpi6->Draw();

  //  c1->Print("kaon_distributions_5_on_50.pdf","pdf");


  //=======================================================================               

  //XBIN 2

  cout << "===== XBIN 2: x=0.15 ========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c2->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //myTree1->Draw("fpi>>myfpix21",totalCut&&xbinCut2);
   myTree1->Draw("fpi>>myfpix21",xbinCut2);

    cout << "fpi mean for Q2=3.75, x=0.15 is: " << myfpix21.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix21->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix21);

  myXpiHist1->GetXaxis()->SetTitle("xkaon");

 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c2->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix22",totalCut&&xbinCut2);
  myTree2->Draw("fpi>>myfpix22",xbinCut2);

  cout << "fpi mean for Q2=12, x=0.15 is: " << myfpix22.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix22->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix22);

  // myfpix22->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c2->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix23",totalCut&&xbinCut2); 
  myTree3->Draw("fpi>>myfpix23",xbinCut2); 

   cout << "fpi mean for Q2=30, x=0.15 is: " << myfpix23.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix23->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix23);

  //myfpix23->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c2->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix24",totalCut&&xbinCut2);
      myTree4->Draw("fpi>>myfpix24",xbinCut2);

    cout << "fpi mean for Q2=60, x=0.15 is: " << myfpix24.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix24->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix24);

  //myfpix24->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c2->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  //myTree5->Draw("fpi>>myfpix25",totalCut&&xbinCut2);
      myTree5->Draw("fpi>>myfpix25",xbinCut2);

    cout << "fpi mean for Q2=120, x=0.15 is: " << myfpix25.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix25->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix25);

  //myfpix25->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c2->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix26",totalCut&&xbinCut2);
   myTree6->Draw("fpi>>myfpix26",xbinCut2);

  cout << "fpi mean for Q2=240, x=0.15 is: " << myfpix26.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix26->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix26);

  //myfpix26->Draw();

  //=======================================================================               

  //XBIN 3

  cout << "====== XBIN 3: x=0.25 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c3->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  // myTree1->Draw("fpi>>myfpix31",totalCut&&xbinCut3);
        myTree1->Draw("fpi>>myfpix31",xbinCut3);

    cout << "fpi mean for Q2=3.75, x=0.25 is: " << myfpix31.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix31->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix31);

  myXpiHist1->GetXaxis()->SetTitle("xkaon");

 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c3->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix32",totalCut&&xbinCut3);
   myTree2->Draw("fpi>>myfpix32",xbinCut3);
  
   cout << "fpi mean for Q2=12, x=0.25 is: " << myfpix32.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix32->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix32);

  // myfpix32->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c3->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix33",totalCut&&xbinCut3); 
  myTree3->Draw("fpi>>myfpix33",xbinCut3); 

   cout << "fpi mean for Q2=30, x=0.25 is: " << myfpix33.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix33->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix33);

  //myfpix33->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c3->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix34",totalCut&&xbinCut3);
       myTree4->Draw("fpi>>myfpix34",xbinCut3);

    cout << "fpi mean for Q2=60, x=0.25 is: " << myfpix34.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix34->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix34);

  //myfpix34->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c3->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  //myTree5->Draw("fpi>>myfpix35",totalCut&&xbinCut3);
         myTree5->Draw("fpi>>myfpix35",xbinCut3);

    cout << "fpi mean for Q2=120, x=0.25 is: " << myfpix35.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix35->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix35);

  //myfpix35->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c3->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix36",totalCut&&xbinCut3);
    myTree6->Draw("fpi>>myfpix36",xbinCut3);

  cout << "fpi mean for Q2=240, x=0.25 is: " << myfpix36.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix36->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix36);

  //myfpix36->Draw();

  //=======================================================================               

  //XBIN 4

  cout << "====== XBIN 4: x=0.35 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c4->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Draw("fpi>>myfpix41",totalCut&&xbinCut4);
         myTree1->Draw("fpi>>myfpix41",xbinCut4);

    cout << "fpi mean for Q2=3.75, x=0.35 is: " << myfpix41.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix41->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix41);


 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c4->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix42",totalCut&&xbinCut4);
   myTree2->Draw("fpi>>myfpix42",xbinCut4);

  cout << "fpi mean for Q2=12, x=0.35 is: " << myfpix42.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix42->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix42);

  // myfpix42->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c4->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  // myTree3->Draw("fpi>>myfpix43",totalCut&&xbinCut4); 
     myTree3->Draw("fpi>>myfpix43",xbinCut4); 

   cout << "fpi mean for Q2=30, x=0.35 is: " << myfpix43.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix43->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix43);

  //myfpix43->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c4->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  //  myTree4->Draw("fpi>>myfpix44",totalCut&&xbinCut4);
      myTree4->Draw("fpi>>myfpix44",xbinCut4);

    cout << "fpi mean for Q2=60, x=0.35 is: " << myfpix44.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix44->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix44);

  //myfpix44->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c4->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  //myTree5->Draw("fpi>>myfpix45",totalCut&&xbinCut4);
    myTree5->Draw("fpi>>myfpix45",xbinCut4);

    cout << "fpi mean for Q2=120, x=0.35 is: " << myfpix45.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix45->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix45);

  //myfpix45->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c4->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix46",totalCut&&xbinCut4);
    myTree6->Draw("fpi>>myfpix46",xbinCut4);

  cout << "fpi mean for Q2=240, x=0.35 is: " << myfpix46.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix46->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix46);

  //myfpix46->Draw();

  //=======================================================================               

  //XBIN 5

  cout << "====== XBIN 5: x=0.45 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c5->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Draw("fpi>>myfpix51",totalCut&&xbinCut5);
        myTree1->Draw("fpi>>myfpix51",xbinCut5);
    cout << "fpi mean for Q2=3.75, x=0.45 is: " << myfpix51.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix51->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix51);


 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c5->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix52",totalCut&&xbinCut5);
    myTree2->Draw("fpi>>myfpix52",xbinCut5);

  cout << "fpi mean for Q2=12, x=0.45 is: " << myfpix52.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix52->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix52);

  // myfpix52->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c5->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix53",totalCut&&xbinCut5); 
    myTree3->Draw("fpi>>myfpix53",xbinCut5); 

   cout << "fpi mean for Q2=30, x=0.45 is: " << myfpix53.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix53->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix53);

  //myfpix53->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c5->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix54",totalCut&&xbinCut5);
      myTree4->Draw("fpi>>myfpix54",xbinCut5);

    cout << "fpi mean for Q2=60, x=0.45 is: " << myfpix54.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix54->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix54);

  //myfpix54->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c5->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  // myTree5->Draw("fpi>>myfpix55",totalCut&&xbinCut5);
        myTree5->Draw("fpi>>myfpix55",xbinCut5);

    cout << "fpi mean for Q2=120, x=0.45 is: " << myfpix55.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix55->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix55);

  //myfpix55->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c5->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix56",totalCut&&xbinCut5);
    myTree6->Draw("fpi>>myfpix56",xbinCut5);

  cout << "fpi mean for Q2=240, x=0.45 is: " << myfpix56.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix56->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix56);

  //myfpix56->Draw();


  //=======================================================================               

  //XBIN 6

  cout << "====== XBIN 6: x=0.55 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c6->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Draw("fpi>>myfpix61",totalCut&&xbinCut6);
        myTree1->Draw("fpi>>myfpix61",xbinCut6);

    cout << "fpi mean for Q2=3.75, x=0.55 is: " << myfpix61.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix61->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix61);


 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c6->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix62",totalCut&&xbinCut6);
    myTree2->Draw("fpi>>myfpix62",xbinCut6);
  cout << "fpi mean for Q2=12, x=0.55 is: " << myfpix62.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix62->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix62);

  // myfpix62->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c6->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix63",totalCut&&xbinCut6); 
    myTree3->Draw("fpi>>myfpix63",xbinCut6); 

   cout << "fpi mean for Q2=30, x=0.55 is: " << myfpix63.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix63->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix63);

  //myfpix63->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c6->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix64",totalCut&&xbinCut6);
      myTree4->Draw("fpi>>myfpix64",xbinCut6);

    cout << "fpi mean for Q2=60, x=0.55 is: " << myfpix64.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix64->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix64);

  //myfpix64->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c6->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  // myTree5->Draw("fpi>>myfpix65",totalCut&&xbinCut6);
        myTree5->Draw("fpi>>myfpix65",xbinCut6);

    cout << "fpi mean for Q2=120, x=0.55 is: " << myfpix65.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix65->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix65);

  //myfpix65->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c6->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix66",totalCut&&xbinCut6);
    myTree6->Draw("fpi>>myfpix66",xbinCut6);

  cout << "fpi mean for Q2=240, x=0.55 is: " << myfpix66.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix66->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix66);

  //myfpix66->Draw();

  //=======================================================================               

  //XBIN 7

  cout << "====== XBIN 7: x=0.65 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c7->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Draw("fpi>>myfpix71",totalCut&&xbinCut7);
        myTree1->Draw("fpi>>myfpix71",xbinCut7);

    cout << "fpi mean for Q2=3.75, x=0.65 is: " << myfpix71.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix71->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix71);


 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c7->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix72",totalCut&&xbinCut7);
    myTree2->Draw("fpi>>myfpix72",xbinCut7);

  cout << "fpi mean for Q2=12, x=0.65 is: " << myfpix72.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix72->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix72);

  // myfpix72->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c7->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix73",totalCut&&xbinCut7); 
    myTree3->Draw("fpi>>myfpix73",xbinCut7); 

   cout << "fpi mean for Q2=30, x=0.65 is: " << myfpix73.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix73->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix73);

  //myfpix73->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c7->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix74",totalCut&&xbinCut7);
      myTree4->Draw("fpi>>myfpix74",xbinCut7);

    cout << "fpi mean for Q2=60, x=0.65 is: " << myfpix74.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix74->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix74);

  //myfpix74->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c7->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  // myTree5->Draw("fpi>>myfpix75",totalCut&&xbinCut7);
        myTree5->Draw("fpi>>myfpix75",xbinCut7);

    cout << "fpi mean for Q2=120, x=0.65 is: " << myfpix75.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix75->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix75);

  //myfpix75->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c7->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix76",totalCut&&xbinCut7);
    myTree6->Draw("fpi>>myfpix76",xbinCut7);

  cout << "fpi mean for Q2=240, x=0.65 is: " << myfpix76.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix76->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix76);

  //myfpix76->Draw();


  //=======================================================================               

  //XBIN 8

  cout << "====== XBIN 8: x=0.75 =========" << endl;

  // XPI FOR 5 ON 50: Q2=3.75 GeV
  c8->cd(1);

  //Full Tree
  TTree* myTree1 = (TTree*) mySimFile1->Get("Evnts");
  //  myTree1->Draw("fpi>>myfpix81",totalCut&&xbinCut8);
        myTree1->Draw("fpi>>myfpix81",xbinCut8);

    cout << "fpi mean for Q2=3.75, x=0.75 is: " << myfpix81.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix81->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix81);


 
  // XPI FOR 5 ON 50: Q2=12 GeV
  c8->cd(2);

  //Full Tree
  TTree* myTree2 = (TTree*) mySimFile2->Get("Evnts");
  //myTree2->Draw("fpi>>myfpix82",totalCut&&xbinCut8);
    myTree2->Draw("fpi>>myfpix82",xbinCut8);

  cout << "fpi mean for Q2=12, x=0.75 is: " << myfpix82.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix82->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix82);

  // myfpix82->Draw();


  // XPI FOR 5 ON 50: Q2=30 GeV
  c8->cd(3);

  //Full Tree
  TTree* myTree3 = (TTree*) mySimFile3->Get("Evnts");
  //myTree3->Draw("fpi>>myfpix83",totalCut&&xbinCut8); 
    myTree3->Draw("fpi>>myfpix83",xbinCut8); 

   cout << "fpi mean for Q2=30, x=0.75 is: " << myfpix83.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix83->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix83);

  //myfpix83->Draw();


  // XPI FOR 5 ON 50: Q2=60 GeV
  c8->cd(4);

  //Full Tree
  TTree* myTree4 = (TTree*) mySimFile4->Get("Evnts");
  // myTree4->Draw("fpi>>myfpix84",totalCut&&xbinCut8);
      myTree4->Draw("fpi>>myfpix84",xbinCut8);

    cout << "fpi mean for Q2=60, x=0.75 is: " << myfpix84.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix84->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix84);

  //myfpix84->Draw();


  // XPI FOR 5 ON 50: Q2=120 GeV
  c8->cd(5);

  //Full Tree
  TTree* myTree5 = (TTree*) mySimFile5->Get("Evnts");
  // myTree5->Draw("fpi>>myfpix85",totalCut&&xbinCut8);
        myTree5->Draw("fpi>>myfpix85",xbinCut8);

    cout << "fpi mean for Q2=120, x=0.75 is: " << myfpix85.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix85->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix85);

  //myfpix85->Draw();


  // XPI FOR 5 ON 50: Q2=240 GeV
  c8->cd(6);

  //Full Tree
  TTree* myTree6 = (TTree*) mySimFile6->Get("Evnts");
  //myTree6->Draw("fpi>>myfpix86",totalCut&&xbinCut8);
    myTree6->Draw("fpi>>myfpix86",xbinCut8);

  cout << "fpi mean for Q2=240, x=0.75 is: " << myfpix86.GetMean() << endl;

  Int_t counts=0;
  counts =  myfpix86->GetEntries();
  cout << "counts from 5 on 50: " << counts << endl;
  fHistoList->AddAtFree(myfpix86);

  //myfpix86->Draw();


  gROOT->Reset();


  return 0;

}
