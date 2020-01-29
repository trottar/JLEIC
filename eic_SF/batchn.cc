void batchn()
{
    gROOT->ProcessLine(".L TDISMC_EICn.cpp+");
    gROOT->ProcessLine("mainx(0.055,0.3,0.1,100,20001)");
}
