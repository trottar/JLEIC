c
      FUNCTION ALPHAsl(Q)
      IMPLICIT REAL*8 (A-H,P-Z)
      pi=2.d0*dasin(1.d0)
      ZETAM = 91.161D0
      ALPHAZ = 0.12D0
      FLAVOR = 4.0D0
      BZERO = 11.0D0 - 2.0D0*FLAVOR/ 3.0D0
      ALPHAsl = ALPHAZ/ (1 + ALPHAZ*BZERO*DLOG(dsqrt(Q)/ZETAM)/2.0D0/PI)
c here use for consistency alphas(SQRT(Q))
c      write(6,*) ' blw alphas ',alphasl
      alphasl = alphas(dsqrt(Q))
c      write(6,*) ' my alphas ',alphasl
      RETURN
      END
