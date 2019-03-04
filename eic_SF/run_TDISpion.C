#include <TProof.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

void run_TDISpion(){

  TChain ch("Evnts");

  ch.Add("TDISpion.root");
  
  TProof *proof = TProof::Open("workers=4");
  //proof->SetProgressDialog(0);  
  ch.SetProof();
  //ch.Process("TDISpion.C+",option);
  ch.Process("TDISpion.C+");
  proof->Close();

}
