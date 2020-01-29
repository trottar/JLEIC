{
  gROOT->Reset();

  c1 = new TCanvas("c1","Projections Kaon",200,10,700,500);

  c1->SetFillColor(10);
  c1->SetGrid();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);

  c1->Divide(3,2);

  Int_t n = 10;
  Float_t x1[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y1[n]  = {1,1,1,1,1,1,1,1,1,1};
  Float_t ex1[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey1[n] = {sqrt(810)/810,sqrt(2105)/2105,sqrt(1226)/1226,sqrt(657)/657,sqrt(448)/448,sqrt(324)/324,sqrt(225)/225,sqrt(184)/184,sqrt(122)/122,sqrt(108)/108};
  gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
  gr1->SetTitle("s=1000 GeV, Q2=3.75 GeV2");
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(21);

  gr1->SetMinimum(0.0);
  gr1->SetMaximum(2.0);
  gr1->GetXaxis()->SetTitle("xkaon");

  c1->cd(1);
  gr1->Draw("AP");


  c1->Update();

  Float_t x2[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y2[n]  = {1,1,1,1,1,1,1,1,1,1};
  Float_t ex2[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey2[n] = {sqrt(1181)/1181,sqrt(2378)/2378,sqrt(1376)/1376,sqrt(748)/748,sqrt(502)/502,sqrt(343)/343,sqrt(270)/270,sqrt(220)/220,sqrt(123)/123,sqrt(138)/138};
  gr2 = new TGraphErrors(n,x2,y2,ex2,ey2);
  gr2->SetTitle("s=1000 GeV, Q2=12 GeV2");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);

  gr2->SetMinimum(0.0);
  gr2->SetMaximum(2.0);
  gr2->GetXaxis()->SetTitle("xkaon");

  c1->cd(2);
  gr2->Draw("AP");

  Float_t x3[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y3[n]  = {1,1,1,1,1,1,1,1,1,1};
  Float_t ex3[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey3[n] = {sqrt(793)/793,sqrt(1811)/1811,sqrt(1168)/1168,sqrt(730)/730,sqrt(495)/495,sqrt(364)/364,sqrt(237)/237,sqrt(182)/182,sqrt(160)/160,sqrt(110)/110};
  gr3 = new TGraphErrors(n,x3,y3,ex3,ey3);
  gr3->SetTitle("s=1000 GeV, Q2=30 GeV2");
  gr3->SetMarkerColor(4);
  gr3->SetMarkerStyle(21);

  gr3->SetMinimum(0.0);
  gr3->SetMaximum(2.0);
  gr3->GetXaxis()->SetTitle("xkaon");

  c1->cd(3);
  gr3->Draw("AP");

  Float_t x4[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y4[n]  = {1,1,1,1,1,1,1,1,1,1};
  Float_t ex4[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey4[n] = {sqrt(463)/463,sqrt(1357)/1357,sqrt(1009)/1009,sqrt(613)/613,sqrt(381)/381,sqrt(225)/225,sqrt(156)/156,sqrt(135)/135,sqrt(80)/80,sqrt(80)/80};
  gr4 = new TGraphErrors(n,x4,y4,ex4,ey4);
  gr4->SetTitle("s=1000 GeV, Q2=60 GeV2");
  gr4->SetMarkerColor(4);
  gr4->SetMarkerStyle(21);

  gr4->SetMinimum(0.0);
  gr4->SetMaximum(2.0);
  gr4->GetXaxis()->SetTitle("xkaon");

  c1->cd(4);
  gr4->Draw("AP");


  Float_t x5[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y5[n]  = {1,1,1,1,1,1,1,1,1,9999};
  Float_t ex5[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey5[n] = {sqrt(1)/1,sqrt(800)/800,sqrt(767)/767,sqrt(481)/481,sqrt(203)/203,sqrt(94)/94,sqrt(39)/39,sqrt(14)/14,sqrt(10)/10,0.};
  gr5 = new TGraphErrors(n,x5,y5,ex5,ey5);
  gr5->SetTitle("s=1000 GeV, Q2=120 GeV2");
  gr5->SetMarkerColor(4);
  gr5->SetMarkerStyle(21);

  gr5->SetMinimum(0.0);
  gr5->SetMaximum(2.0);
  gr5->GetXaxis()->SetTitle("xkaon");

  c1->cd(5);
  //  gr5->Draw("ALP"); - IF WANT TO DRAW WITH LINE BETWEEN POINTS
  gr5->Draw("AP");

  Float_t x6[n]  = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75,0.85,0.95};
  Float_t y6[n]  = {9999,9999,1,1,1,9999,9999,9999,9999,9999};
  Float_t ex6[n] = {.005,.001,.007,.007,.004,.005,.006,.007,.008,.005};
  Float_t ey6[n] = {sqrt(0)/0,sqrt(0)/0,sqrt(328)/328,sqrt(215)/215,sqrt(10)/10,sqrt(0)/0,sqrt(0)/0,sqrt(0)/0,sqrt(0)/0,sqrt(0)/0};
  gr6 = new TGraphErrors(n,x6,y6,ex6,ey6);
  gr6->SetTitle("s=1000 GeV, Q2=240 GeV2");
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);

  gr6->SetMinimum(0.0);
  gr6->SetMaximum(2.0);
  gr6->GetXaxis()->SetTitle("xkaon");

  c1->cd(6);
  gr6->Draw("AP");

  c1->Print("kaon_projections_5_on_50.pdf","pdf");

}
