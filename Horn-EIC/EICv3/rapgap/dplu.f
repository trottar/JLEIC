*CMZ :  2.07/03 15/05/99  14.52.10  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
C*********************************************************************

      FUNCTION DPLU(I,J)
      IMPLICIT NONE
C...Purpose: to provide various real-valued event related data.
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
*KEEP,RGLUPARM.
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


*KEND.
      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      REAL ULMASS
	Double precision DLANGL,PSUM,DPLU,PR,PMR
	Integer I,J,I1,J1,LUCHGE
      EXTERNAL DLANGL,ULMASS,LUCHGE
C      SAVE

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      SAVE /LUJETS/
      DIMENSION PSUM(4)
C...Set default value. For I = 0 sum of momenta or charges,
C...or invariant mass of system.
      DPLU=0.D0
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN
         DO 10 I1=1,N
   10    IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) DPLU=DPLU+P(I1,J)
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN
         DO 20 J1=1,4
            PSUM(J1)=0.D0
            DO 20 I1=1,N
   20    IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+P(I1,J1)
         DPLU=DSQRT(DMAX1(0.D0, PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-
     +   PSUM(3)**2))
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN
         DO 30 I1=1,N
   30    IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) DPLU=DPLU+
     +   DBLE(LUCHGE(K(I1,2)))/3.D0
      ELSEIF(I.EQ.0) THEN

C...Direct readout of P matrix.
      ELSEIF(J.LE.5) THEN
         DPLU=P(I,J)

C...Charge, total momentum, transverse momentum, transverse mass.
      ELSEIF(J.LE.12) THEN
         IF(J.EQ.6) DPLU=DBLE(LUCHGE(K(I,2)))/3.D0
         IF(J.EQ.7.OR.J.EQ.8) DPLU=P(I,1)**2+P(I,2)**2+P(I,3)**2
         IF(J.EQ.9.OR.J.EQ.10) DPLU=P(I,1)**2+P(I,2)**2
         IF(J.EQ.11.OR.J.EQ.12) DPLU=P(I,5)**2+P(I,1)**2+P(I,2)**2
         IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) DPLU=DSQRT(DPLU)

C...Theta and phi angle in radians or degrees.
      ELSEIF(J.LE.16) THEN
         IF(J.LE.14) DPLU=DLANGL(P(I,3),DSQRT(P(I,1)**2+P(I,2)**2))
         IF(J.GE.15) DPLU=DLANGL(P(I,1),P(I,2))
         IF(J.EQ.14.OR.J.EQ.16) DPLU=DPLU*DBLE(180./PARU(1))

C...True rapidity, rapidity with pion mass, pseudorapidity.
      ELSEIF(J.LE.19) THEN
         PMR=0.D0
         IF(J.EQ.17) PMR=P(I,5)
         IF(J.EQ.18) PMR=DBLE(ULMASS(211))
         PR=DMAX1(1D-20,PMR**2+P(I,1)**2+P(I,2)**2)
         DPLU=SIGN(DLOG(DMIN1((DSQRT(PR+(P(I,3)**2)) +ABS(P(I,3))
     +   )/DSQRT(PR),1D20)),P(I,3))

C...Energy and momentum fractions (only to be used in CM frame).
      ELSEIF(J.LE.25) THEN
         IF(J.EQ.20) DPLU=2.D0*DSQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)/
     +   DBLE(PARU(21))
         IF(J.EQ.21) DPLU=2.D0*P(I,3)/DBLE(PARU(21))
         IF(J.EQ.22) DPLU=2.D0*DSQRT(P(I,1)**2+P(I,2)**2)/DBLE(PARU(21))
         IF(J.EQ.23) DPLU=2.D0*P(I,4)/DBLE(PARU(21))
         IF(J.EQ.24) DPLU=(P(I,4)+P(I,3))/DBLE(PARU(21))
         IF(J.EQ.25) DPLU=(P(I,4)-P(I,3))/DBLE(PARU(21))
      ENDIF

      RETURN
      END
