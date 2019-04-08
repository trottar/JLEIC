void batch()
{
    gROOT->ProcessLine(".L TDISMC_EIC.cpp+");
    gROOT->ProcessLine("mainx(0.0001,1.,1.0,50.0,101)");
}
