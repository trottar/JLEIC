*CMZ :  2.06/39 16/07/98  11.40.09  by  Hannes Jung
*CMZ :  2.06/38 10/07/98  19.31.42  by  Hannes Jung
*CMZ :  2.06/20 26/11/97  09.43.14  by  Hannes Jung
*CMZ :  2.06/19 19/11/97  10.41.57  by  Hannes Jung
*CMZ :  2.06/17 12/08/97  20.27.32  by  Julian Phillips
*-- Author :
*-- Author :    Julian Phillips   05/08/97
*###########################################################
* Main routine for parameterisations of final H1 1994
* QCD fits to F2D(3)
*
*   Input:    zt=x_{i/IP}
*            q2t=photon virtuality
*           ifit=1,2,3 (NLO fits as in paper)
*           ifit=4,5,6 (LO fits)
*           Call with IFIT>0 to initialise the parameterisation,
*           then call with ifit=0 for subsequent access.
*
*   Output: xpq(-6:6): PDG style array of partons, 0=gluon
*
* Note that the parameterisations are leading order, and are
* only valid in the kinematic range
*         0.04<beta<1.0
*         4.5<Q2<75 GeV2
* THIS IS THE KINEMATIC RANGE IN WHICH RESULTS ARE AVAILABLE
* TO THE WORLD.
*###########################################################
* To aid H1 Monte Carlo production, I have extended the
* parameterisation to the largest kinematic range in
* which I have any faith whatsoever in the parameterisation:
*         0.01<beta<1.0
*         3.0<Q2<100 GeV2
* For values of Q^2 or beta outside the specified range,
* the parameterisation is taken from the nearest valid
* beta and/or q2 value.
*
* YOU SHOULD NOT EXPECT THE MONTE CARLO TO LOOK LIKE THE
* DATA IN THE REGION Q2<4.5 OR BETA<0.04 SO DON'T SAY
* I DIDN'T WARN YOU.
*
* The pion structure function used is just an unmodified
* one, but RAPGAP will generate events with as if the
* particle was neutral, but with more u and d quarks
* in the current jet than expected.
*
* J.P.Phillips, 5/8/97
*###########################################################
      SUBROUTINE QCD_1994(ZT1,Q2T1,XPQ,IFIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IQ2,IZ,NQ2,NZ,IFIT
      PARAMETER (NQ2MAX=30,NZMAX=50)
      DIMENSION IPAR(3), XPQ(-6:6)
      DATA IPAR/1,4,0/
      DIMENSION UDS_GRD(NZMAX,NQ2MAX),C_GRD(NZMAX,NQ2MAX),
     +          G_GRD(NZMAX,NQ2MAX)
      COMMON /CPARAM/ UDS_GRD,C_GRD,G_GRD

C Initialise Parameterisation from data statements
      IF(IFIT.GT.0) THEN
         WRITE(6,*)'Initialising fit',IFIT
         IF(IFIT.EQ.1) CALL I_NLO_Q3G0
         IF(IFIT.EQ.2) CALL I_NLO_Q3G1
         IF(IFIT.EQ.3) CALL I_NLO_Q3G3
         IF(IFIT.EQ.4) CALL I_LO_Q3G0
         IF(IFIT.EQ.5) CALL I_LO_Q3G1
         IF(IFIT.EQ.6) CALL I_LO_Q3G3
      ENDIF

C Initialise xpq
      DO I=-6,6
         XPQ(I)=0.D0
      ENDDO
      ZT = ZT1
      Q2T = Q2T1
      IF(ZT.GT.0.99D0) THEN
c            WRITE(6,*)' QCD_94 zt,q2t ',ZT,Q2T
         ZT = 0.99D0
      ENDIF
      IF(Q2T.GT.99.9D0) THEN
c            WRITE(6,*)' QCD_94 zt,q2t ',ZT,Q2T
         Q2T = 99.9D0
      ENDIF

C Define grid
*! Q2min
      Q2L = 3.D0
*! Q2max
      Q2U = 99.99D0
*! beta_min
      ZL  = 0.01D0
*! beta_mid (end of log binning)
      ZM  = 0.50D0
*! beta_max
      ZU  = 0.9999D0
*! Points in Q2 grid
      NQ2 = 30
*! Points in log beta grid
      NZa = 25
*! Total points in grid
      NZ  = 50
      DQ2 = DLOG(Q2U/Q2L)/DFLOAT(NQ2-1)
      DZa = DLOG(ZM/ZL)/DFLOAT(NZa-1)
      DZb = (ZU-ZM)/DFLOAT(NZ-NZa)

C Check input parameters
      Q2T=MIN(Q2U-0.01D0,MAX(Q2L,Q2T))
      ZT =MIN(ZU,MAX(ZL,ZT))

C Lower Grid Point
      IQ2 = INT(DLOG(Q2T/Q2L)/DQ2)+1
      IZ  = INT(DLOG(MIN(ZT,ZM)/ZL)/DZa)+1
     +      + INT(MAX(0.D0,ZT-ZM)/DZb)

C Central and boundary z,Q2 points for interpolation
      DLQ2  = DLOG(Q2T)
      DLQ21 = DLOG(Q2L)+DFLOAT(IQ2-1)*DQ2
      DLQ22 = DLOG(Q2L)+DFLOAT(IQ2)*DQ2
      IF(IZ.LT.NZa) THEN
*! Logarithmic interpolation
         DLZ = DLOG(ZT)
         DLZ1 = DLOG(DEXP(DLOG(ZL)+DFLOAT(MIN(IZ,NZa)-1)*DZa) +
     +   DFLOAT(MAX(IZ-NZa,0))*DZb)
         DLZ2 = DLOG(DEXP(DLOG(ZL)+DFLOAT(MIN(IZ+1,NZa)-1)*DZa) +
     +   DFLOAT(MAX(IZ+1-NZa,0))*DZb)
      ELSE
*! Linear interpolation
         DLZ = ZT
         DLZ1 = DEXP(DLOG(ZL)+DFLOAT(MIN(IZ,NZa)-1)*DZa) +
     +   DFLOAT(MAX(IZ-NZa,0))*DZb
         DLZ2 = DEXP(DLOG(ZL)+DFLOAT(MIN(IZ+1,NZa)-1)*DZa) +
     +   DFLOAT(MAX(IZ+1-NZa,0))*DZb
      ENDIF

C Light Flavour Singlet
      UDS11 = UDS_GRD(IZ,IQ2)
      UDS12 = UDS_GRD(IZ,IQ2+1)
      UDS21 = UDS_GRD(IZ+1,IQ2)
      UDS22 = UDS_GRD(IZ+1,IQ2+1)

C Charm
      C11 = C_GRD(IZ,IQ2)
      C12 = C_GRD(IZ,IQ2+1)
      C21 = C_GRD(IZ+1,IQ2)
      C22 = C_GRD(IZ+1,IQ2+1)

C Gluon
      G11 = G_GRD(IZ,IQ2)
      G12 = G_GRD(IZ,IQ2+1)
      G21 = G_GRD(IZ+1,IQ2)
      G22 = G_GRD(IZ+1,IQ2+1)

C Interpolate to parton distributions to ZT,Q2T
      UDS_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     +     UDS11,UDS12,UDS21,UDS22)
      C_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     +     C11,C12,C21,C22)
      G_T = XYINTER(DLZ1,DLZ2,DLQ21,DLQ22,DLZ,DLQ2,
     +     G11,G12,G21,G22)

C Fill XPQ array
      XPQ(-4) = C_T
      XPQ(-3) = UDS_T
      XPQ(-2) = UDS_T
      XPQ(-1) = UDS_T
      XPQ( 0) = G_T
      XPQ( 1) = UDS_T
      XPQ( 2) = UDS_T
      XPQ( 3) = UDS_T
      XPQ( 4) = C_T

C Check output
      DO I=-6,6
         IF(XPQ(I).LT.-0.1D-10) THEN
            WRITE(6,*)'Error in H1 fit parameterisation',XPQ(I)
            WRITE(6,*)' at zt,q2t ',ZT,Q2T
            STOP
         ENDIF
      ENDDO

      RETURN
      END
