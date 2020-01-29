#! /usr/bin/perl
use Math::Trig;
use Math::Complex;

my ($ir, $iqq, $ixx) = @ARGV;       

# This kinematics is what we decide as final kinematic bins for pseudo-data generation @ Aug.10,2016
# this is for the TDIS(Tagged Deeply Inelastic Scattering) MC simulation with CTEQ DIS cross-section
# parametriztion in the frame work of Jlab Electron0-Ion-Collder
#
# This code is a sort of extension of STEG but there is significant modification for TDIS
#
                                           
$In_xMin=0.1;
$In_xMax=0.3;
$In_Q2Min=1.;
$In_Q2Max=100.;

$NumberofBinXbj = 1;
$NumberofBinQ2 = 1;
$dlogxBj = ($In_xMax)-log($In_xMin);
$ddlogxBj = abs($dlogxBj)/$NumberofBinXbj;
$dlogQ2 = ($In_Q2Max)-($In_Q2Min);
$ddlogQ2 = abs($dlogQ2)/$NumberofBinQ2;

$srunnum = 10000000+100000*($iqq+1)+1000*$ixx;
$erunnum = 10000000+100000*($iqq+1)+1000*$ixx+1;
                                

$q1 = ($In_Q2Min)+$ddlogQ2*$iqq;
$q2 = ($In_Q2Min)+$ddlogQ2*($iqq+1);
                         
$x1 = ($In_xMin)+$ddlogxBj*$ixx;
$x2 = ($In_xMin)+$ddlogxBj*($ixx+1);

$ir = $srunnum;

system("rm -f TDISMC_EIC_cpp.d TDISMC_EIC_cpp.so");
system("rm -f TDIS_lund.txt");

system("rm -f batch.cc");
# Create jbatch
open(file_handle,">batch.cc");
print file_handle <<END_OF_FILE1;
void batch()
{
    gROOT->ProcessLine(".L TDISMC_EIC.cpp+");
    gROOT->ProcessLine("mainx($x1,$x2,$q1,$q2,$ir)");
}
END_OF_FILE1
    close(file_handle);

system("root -b -q batch.cc");
#system("scp TDIS_lund.txt ./jleic_step1.temp");
print "currently reading : RunNumber: $ir,     xBJ= $x1, $x2,    Q2= $q1,  $q2 \n";
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ run number $ir is completed......\n\n\n";
#system("ls -lt jleic_step1.temp");


# This is the end of script.
