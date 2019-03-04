void batch()
{
    gROOT->ProcessLine(".L TDISMC_EIC.cpp+");
    gROOT->ProcessLine("mainx(0.001,1.,1,100,101)");
}
