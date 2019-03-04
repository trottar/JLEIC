#!/bin/csh -f

@ i = 1
@ n = 20
while ( $i <= $n )
  echo Starting $i
  setenv RGFIXED_PE -6.0
  set outname=rgfixed-060-$i
  setenv RGFIXED_MCSEED mcseed.dat
  setenv RGFIXED_OUTPUT $outname

  ./rgfixed.x | grep -i sigma >& $outname.log

  @ i = $i + 1

end
