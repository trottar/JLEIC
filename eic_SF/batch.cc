void batch()
{
    gROOT->ProcessLine(".L TDISMC_EIC.cpp+");
    gROOT->ProcessLine("mainx(0.00001,1.,1.0,300.0,103)");
}
