#!/bin/csh -f

@ i = 1
#@ n = 20
@ n = 1
while ( $i <= $n )
  echo Starting $i
  setenv RGFIXED_PE -6.0
  set outname=rgfixed-060-$i
  setenv RGFIXED_MCSEED mcseed.dat
  setenv RGFIXED_OUTPUT $outname

#  ./rgfixed.x | grep -i sigma >& $outname.log
  ./rgfixed.x >& $outname.log
  grep -i sigma $outname.log >& $outname.sigma


  @ i = $i + 1

end
