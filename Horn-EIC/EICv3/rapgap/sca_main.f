*CMZ :  2.08/05 28/03/2000  14.07.14  by  Frank-Peter Schilling
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*-- Author :    Julian Phillips   13/09/96

*     Semiclassical Model Main Routine

      SUBROUTINE SCA_MAIN(BETA_IN,Q2_IN,XPQ,X_POM_IN,T2_IN)

      IMPLICIT REAL*8 (A-G,O-Z)
      REAL*4 T2_IN,X_POM_IN,BETA_IN,Q2_IN,XPQ(-6:6)
      DOUBLE PRECISION XPQD(-6:6)

      INTEGER ICALL
      DATA ICALL/0/
      SAVE ICALL

C     Inform user of folly

      IF(ICALL.EQ.0) THEN
         ICALL=1
         WRITE(6,*)'#############################################'
         WRITE(6,*)'# Semiclassical Model Selected              #'
         WRITE(6,*)'#############################################'
      ENDIF

C     Input quantities are REAL*4 -> Convert to REAL*8

      BETA =DBLE(BETA_IN)
      if(BETA.LT.0.01) BETA = 0.01D0
      if(BETA.GT.0.99) BETA = 0.99D0

      Q2   =DBLE(Q2_IN)
      if(Q2.GT.200.0d0) Q2=200.0d0
      if(Q2.LT.2.0d0)   Q2=2.0d0

      X_POM=DBLE(X_POM_IN)
      T2   =DBLE(T2_IN)

c     Get diffractive PDF's for xpom=0.003

      CALL SCA_PDF(BETA,Q2,XPQD)

c     transform to from xpom=0.003 to input xpom
c     add exponential t-dependence


      parl = 8.16d0
      xpom0 = 0.003d0
      b0=6.0d0

      xpomfac = ((parl-dlog(x_pom))**2/(parl-dlog(xpom0))**2)/(x_pom)
      tfac = DEXP(-(B0*DABS(T2)))*b0

      DO i=-6,6
         xpqd(i)=xpqd(i)*xpomfac*tfac
      enddo

c     convert to real*4

      DO I=-6,6
         XPQ(I)=0.
         IF(ABS(I).LT.5) THEN
           XPQ(I)=REAL(XPQD(I))
         ENDIF
         IF(XPQ(I).LT.1E-10) XPQ(I)=0.
      ENDDO

c     that's it

      RETURN
      END
