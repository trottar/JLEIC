*CMZ :  2.08/02 04/11/99  16.39.08  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  15.04.23  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
C*********************************************************************

      SUBROUTINE DUROBO(DTHE,DPHI,BEX,BEY,BEZ)

C...Purpose: to perform rotations and boosts.
      IMPLICIT NONE
      REAL SP,V
      DOUBLE PRECISION P
      Integer IMIN,IMAX,IMI,IMA,I,J
      Integer N,K
*KEEP,RGLUPARM.
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


*KEND.
      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      SAVE /LUJETS/
      double precision DROT,DPR,VR,DP,DV,DTHE,DPHI,BEX,BEY,BEZ
      double precision DBEX,DBEY,DBEZ
      double precision DBX,DBY,DBZ,DB,DGA,DBP,DGABP,DBV,DGABV
      DIMENSION DROT(3,3),DPR(3),VR(3),DP(4),DV(4)
      REAL SNGL
C...Find range of rotation/boost. Convert boost to double precision.
      IMIN=1
      IF(MSTU(1).GT.0) IMIN=MSTU(1)
      IMAX=N
      IF(MSTU(2).GT.0) IMAX=MSTU(2)
      DBX=BEX
      DBY=BEY
      DBZ=BEZ
      GOTO 20

C...Entry for specific range and double precision boost.
      ENTRY DUDBRB(IMI,IMA,DTHE,DPHI,DBEX,DBEY,DBEZ)
      IMIN=IMI
      IF(IMIN.LE.0) IMIN=1
      IMAX=IMA
      IF(IMAX.LE.0) IMAX=N
      DBX=DBEX
      DBY=DBEY
      DBZ=DBEZ

C...Optional resetting of V (when not set before.)
      IF(MSTU(33).NE.0) THEN
         DO 10 I=MIN(IMIN,MSTU(4)),MIN(IMAX,MSTU(4))
            DO 10 J=1,5
   10    V(I,J)=0.
         MSTU(33)=0
      ENDIF

C...Check range of rotation/boost.
   20 IF(IMIN.GT.MSTU(4).OR.IMAX.GT.MSTU(4)) THEN
         CALL LUERRM(11,'(DUROBO:) range outside LUJETS memory')
         RETURN
      ENDIF

C...Rotate, typically from z axis to direction (theta,phi).
      IF(DTHE**2+DPHI**2.GT.1D-20) THEN
         DROT(1,1)=DCOS(DTHE)*DCOS(DPHI)
         DROT(1,2)=-DSIN(DPHI)
         DROT(1,3)=DSIN(DTHE)*DCOS(DPHI)
         DROT(2,1)=DCOS(DTHE)*DSIN(DPHI)
         DROT(2,2)=DCOS(DPHI)
         DROT(2,3)=DSIN(DTHE)*DSIN(DPHI)
         DROT(3,1)=-DSIN(DTHE)
         DROT(3,2)=0.D0
         DROT(3,3)=DCOS(DTHE)
         DO 50 I=IMIN,IMAX
            IF(K(I,1).LE.0) GOTO 50
            DO 30 J=1,3
               DPR(J)=P(I,J)
   30       VR(J)=DBLE(V(I,J))
            DO 40 J=1,3
               P(I,J)=DROT(J,1)*DPR(1)+DROT(J,2)*DPR(2)+DROT(J,3)*
     +         DPR(3)
   40       V(I,J)=
     +      SNGL(DROT(J,1)*VR(1)+DROT(J,2)*VR(2)+DROT(J,3)*VR(3))
   50    CONTINUE
      ENDIF

C...Boost, typically from rest to momentum/energy=beta.
      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
         DB=DSQRT(DBX**2+DBY**2+DBZ**2)
         IF(DB.GT.0.99999999D0) THEN
C...Rescale boost vector if too close to unity.
            CALL LUERRM(3,'(DUROBO:) boost vector too large')
            call dulist(1)
            DBX=DBX*(0.99999999D0/DB)
            DBY=DBY*(0.99999999D0/DB)
            DBZ=DBZ*(0.99999999D0/DB)
            DB=0.99999999D0
         ENDIF
         DGA=1D0/DSQRT(1D0-DB**2)
         DO 70 I=IMIN,IMAX
            IF(K(I,1).LE.0) GOTO 70
            DO 60 J=1,4
               DP(J)=P(I,J)
   60       DV(J)=DBLE(V(I,J))
            DBP=DBX*DP(1)+DBY*DP(2)+DBZ*DP(3)
            DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP(4))
            P(I,1)=DP(1)+DGABP*DBX
            P(I,2)=DP(2)+DGABP*DBY
            P(I,3)=DP(3)+DGABP*DBZ
            P(I,4)=DGA*(DP(4)+DBP)
            DBV=DBX*DV(1)+DBY*DV(2)+DBZ*DV(3)
            DGABV=DGA*(DGA*DBV/(1D0+DGA)+DV(4))
            V(I,1)=SNGL(DV(1)+DGABV*DBX)
            V(I,2)=SNGL(DV(2)+DGABV*DBY)
            V(I,3)=SNGL(DV(3)+DGABV*DBZ)
            V(I,4)=SNGL(DGA*(DV(4)+DBV))
   70    CONTINUE
      ENDIF

      RETURN
      END
