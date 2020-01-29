void batchk()
{
    gROOT->ProcessLine(".L TDISMC_EICK.cpp+");
    gROOT->ProcessLine("mainxk(0.001,1.,1,100,0)");
}
