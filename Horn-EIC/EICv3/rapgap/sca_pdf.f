*CMZ :  2.08/05 28/03/2000  14.07.14  by  Frank-Peter Schilling
*-- Author :

      SUBROUTINE SCA_PDF(ZT1,Q2T1,XPQ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IQ2,IZ,NQ2,NZ,IFIT

      INTEGER ICALL
      DATA ICALL /0/
      SAVE ICALL

      PARAMETER (NQ2MAX=30,NZMAX=99)

      DIMENSION XPQ(-6:6)

      DIMENSION UDS_GRD(NZMAX,NQ2MAX),C_GRD(NZMAX,NQ2MAX),
     +          G_GRD(NZMAX,NQ2MAX)

      COMMON /CPARAM2/ UDS_GRD,G_GRD,C_GRD

C     Initialise Parameterisation from data statements

      IF(ICALL.EQ.0) THEN
         ICALL=ICALL+1
         WRITE(6,*)' Initialising Semiclassical Model'

         CALL SCA_DAT

      ENDIF

C     Initialise xpq

      DO I=-6,6
         XPQ(I)=0.D0
      ENDDO
      ZT = ZT1
      Q2T = Q2T1
      IF(ZT.GT.0.99D0) THEN
         ZT = 0.99D0
      ENDIF
      IF(Q2T.GT.200.0D0) THEN
         Q2T = 200.D0
      ENDIF

c     define grid

      Q2L = 2.D0       ! q2 min
      Q2U = 200.0D0    ! q2 max
      ZL  = 0.01D0     ! beta min
      ZU  = 0.99D0     ! beta max
      NQ2 = 30         ! Points in Q2 grid
      NZ  = 99         ! Points in beta grid

c     binwidth in q2 and beta

      DQ2 = DLOG(Q2U/Q2L)/DFLOAT(NQ2-1)
      DZ = (ZU-ZL)/DFLOAT(NZ-1)

C     Check input parameters

      Q2T=MIN(Q2U-0.01D0,MAX(Q2L,Q2T))
      ZT =MIN(ZU,MAX(ZL,ZT))

C     Lower Grid Point
      IQ2 = INT(DLOG(Q2T/Q2L)/DQ2)+1
      IZ  = INT(MAX(0.D0,ZT-ZL)/DZ)+1

C     Central and boundary z,Q2 points for interpolation

      DLQ2  = DLOG(Q2T)
      DLQ21 = DLOG(Q2L)+DFLOAT(IQ2-1)*DQ2
      DLQ22 = DLOG(Q2L)+DFLOAT(IQ2)*DQ2

      DLZ = ZT
      DLZ1 = DFLOAT(MAX(IZ,0))*DZ
      DLZ2 = DFLOAT(MAX(IZ+1,0))*DZ

C     Light Flavour Singlet

      UDS11 = UDS_GRD(IZ,IQ2)
      UDS12 = UDS_GRD(IZ,IQ2+1)
      UDS21 = UDS_GRD(IZ+1,IQ2)
      UDS22 = UDS_GRD(IZ+1,IQ2+1)

C     Charm

      C11 = C_GRD(IZ,IQ2)
      C12 = C_GRD(IZ,IQ2+1)
      C21 = C_GRD(IZ+1,IQ2)
      C22 = C_GRD(IZ+1,IQ2+1)

C     Gluon

      G11 = G_GRD(IZ,IQ2)
      G12 = G_GRD(IZ,IQ2+1)
      G21 = G_GRD(IZ+1,IQ2)
      G22 = G_GRD(IZ+1,IQ2+1)

C     Interpolate to parton distributions to ZT,Q2T

      UDS_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     &        UDS11,UDS12,UDS21,UDS22)
      C_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     &      C11,C12,C21,C22)
      G_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     &      G11,G12,G21,G22)

C    Fill XPQ array

      XPQ(-4) = C_T
      XPQ(-3) = UDS_T/6.0d0
      XPQ(-2) = UDS_T/6.0d0
      XPQ(-1) = UDS_T/6.0d0
      XPQ( 0) = G_T
      XPQ( 1) = UDS_T/6.0d0
      XPQ( 2) = UDS_T/6.0d0
      XPQ( 3) = UDS_T/6.0d0
      XPQ( 4) = C_T

C    Check output

      DO I=-6,6
         IF(XPQ(I).LT.-0.1D-10) THEN
            WRITE(6,*)'Error in SCA fit parameterisation',XPQ(I)
            WRITE(6,*)' at zt,q2t ',ZT,Q2T
            STOP
         ENDIF
      ENDDO

c     that's it

      RETURN
      END
