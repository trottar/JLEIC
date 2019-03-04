*KEEP,RGLUJETS.
      include 'rgluparm.inc'

      INTEGER N,K(lupan,5)
      REAL SP(lupan,5),V(lupan,5)
      real*8 P(lupan,5)


      COMMON /LUJETS/ N, K, SP, V
      COMMON /DUJETS/ P

      REAL ULMASS
      real*8 DOT1,DPLU,DLANGL
      EXTERNAL DLANGL,ULMASS,DPLU,LUCHGE,DOT1

*KEND.
