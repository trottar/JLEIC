void batch()
{
    gROOT->ProcessLine(".L TDISMC_EIC.cpp+");
    // gROOT->ProcessLine("mainx(0.00001,1.,1.0,10000.0,103)");
    gROOT->ProcessLine("mainx(0.001,0.008,1.0,8.0,103)");
    // gROOT->ProcessLine("mainx(0.01,0.08,10.0,80.0,103)");
    // gROOT->ProcessLine("mainx(0.1,0.8,100.0,800.0,103)");
}
