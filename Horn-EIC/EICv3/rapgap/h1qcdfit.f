*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 11/02/99  18.53.12  by  Hannes Jung
*CMZ :  2.06/41 21/07/98  16.43.12  by  Hannes Jung
*CMZ :  2.06/34 18/05/98  22.50.06  by  Hannes Jung
*CMZ :  2.06/32 07/04/98  08.28.39  by  Hannes Jung
*CMZ :  2.06/29 15/03/98  13.45.09  by  Hannes Jung
*CMZ :  2.06/28 10/03/98  17.08.43  by  Hannes Jung
*CMZ :  2.06/24 14/01/98  19.40.49  by  Hannes Jung
*CMZ :  2.06/23 29/12/97  17.48.30  by  Hannes Jung
*CMZ :  2.06/18 07/11/97  15.57.13  by  Hannes Jung
*CMZ :  2.06/17 03/11/97  13.25.29  by  Julian Phillips
*CMZ :  2.05/25 04/07/97  12.23.37  by  Hannes Jung
*CMZ :  2.05/20 16/06/97  16.12.29  by  Hannes Jung
*CMZ :  2.03/12 11/10/96  13.36.58  by  Hannes Jung
*CMZ :  2.03/09 04/10/96  16.10.04  by  Julian Phillips
*CMZ :  2.03/08 02/10/96  21.26.02  by  Julian Phillips
*-- Author :    Julian Phillips   13/09/96
*###################################################################
* H1 QCD/Regge Fits
*
* Note on xpq: elements 0,1,2,3,4,5,6 are g,d,u,s,c,b,t
*              ie xpq(1)=d !
* Options
* -------
*
* NPOM = -10 H1 FINAL 1994 FIT (Flux Part)
* NG   = -10 H1 FINAL 1994 FIT (Fit 1, Quarks Only)
* NG   = -11 H1 FINAL 1994 FIT (Fit 2, Quarks + Flat Gluon)
* NG   = -12 H1 FINAL 1994 FIT (Fit 3, Quarks + Flat Gluon)
*####################################################################

      SUBROUTINE H1QCDFIT(BETA_IN,Q2_IN,XPQ,X_POM_IN,T2_IN)

      IMPLICIT REAL*8 (A-G,O-Z)
      REAL*4 T2_IN,X_POM_IN,BETA_IN,Q2_IN,XPQ(-6:6)
      DIMENSION XPQP(-6:6),XPQM(-6:6)
      LOGICAL FIRST
      INTEGER ICALL
      DATA ICALL/0/
      LOGICAL FLAG
      DATA FLAG/.FALSE./
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.

C Function for integral over t
      T_INT(A_TMIN,A_TMAX,B)=(-DEXP(-B*A_TMAX)/B) + (DEXP(-B*A_TMIN)/B)
      DATA FIRST/.TRUE./
      Icall = icall + 1
      if(icall.eq.1) then
         first=.true.
      else
         first=.false.
      endif
C      IF(FIRST) WRITE(6,*)'NPOM',NPOM
C Set Parameters for fits
      IF(NPOM.GT.-10) THEN
*! alpha_prime for the pomeron
         APPOM=0.0D0
*! b0 for pomeron
         B0POM=6.0D0
*! alpha_prime for the meson
         APMES=0.86D0
*! b0 for meson
         B0MES=1.4D0
      ELSEIF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         APPOM=0.26D0
         B0POM=4.6D0
         APMES=0.9D0
         B0MES=2.0D0
      ELSE
         WRITE(6,*)'H1QCDFIT: Unknown fit',NPOM
         STOP
      ENDIF

      IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         CP=1.D0
*! Note: 2*alpha(0)-1 -> this is alpha(0)=1.20
         PP=1.40D0
*!  "        "             "   " alpha(0)=0.57
         PM=0.14D0
         IF(NG.EQ.-10) THEN
            CM=0.011816D0
         ELSEIF(NG.EQ.-11) THEN
            CM=0.0086814D0
         ELSEIF(NG.EQ.-12) THEN
            CM=0.00863D0
         ELSEIF(NG.EQ.-13) THEN
            CM=0.011816D0
         ELSEIF(NG.EQ.-14) THEN
            CM=0.0086814D0
         ELSEIF(NG.EQ.-15) THEN
            CM=0.00863D0
         ELSE
            WRITE(6,*)'H1QCDFIT - Unknown 1994 Final fit',ABS(NG)-9
            STOP
         ENDIF
      ENDIF

C Inform user of folly
      IF(FIRST) THEN
         WRITE(6,*)'#############################################'
         WRITE(6,*)'#           H1QCD fit Selected              #'
         WRITE(6,*)'#       for Q2 > 75  fix pdf at Q2 = 75     #'
         WRITE(6,*)'#       for Q2 < 4.5 fix pdf at Q2 = 4.5    #'
         WRITE(6,*)'#   for beta < 0.04  fix pdf at beta = 0.04 #'
         IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
            WRITE(6,*)'# Final 1994 Fits                           #'
            WRITE(6,*)'# FIT number in paper is ',ABS(NG)-9
            IF(NPOM.EQ.-10) WRITE(6,*)'# Pomeron Part Only'
            IF(NPOM.EQ.-11) WRITE(6,*)'# Meson Part Only'
            IF(NPOM.EQ.-12) WRITE(6,*)'# Pomeron + Meson'
         ELSE
            WRITE(6,*)'# NPOM=',NPOM,'not implemented              #'
            STOP
         ENDIF
         WRITE(6,*)'#############################################'
      ENDIF

C Input quantities are REAL*4 -> Convert to REAL*8
      BETA =DBLE(BETA_IN)
      if(BETA.LT.0.04) BETA = 0.04D0
      Q2   =DBLE(Q2_IN)
      if(Q2.GT.75d0) Q2=75.d0
      if(Q2.LT.4.5d0) Q2=4.5d0
      X_POM=DBLE(X_POM_IN)
      T2   =DBLE(T2_IN)
      IF(NPOM.NE.-3) THEN
C Initialise PDF Pion
         IF(NPOM.EQ.-10 .OR. NPOM.EQ.-11 .OR. NPOM.EQ.-12) THEN
*! LO-Owens Set 1
            CALL PION_INIT(2.D0,1.D0,1.D0)
         ENDIF
      ENDIF

C Final 1994 Fits from paper
      IF(NPOM.LE.-10 .AND. NPOM.GE.-12) THEN
         IF(FIRST) THEN
            IFIT=ABS(NG)-9
         ELSE
            IFIT=0
         ENDIF
         CALL QCD_1994(BETA,Q2,XPQP,IFIT)
c use Q instead of Q2 for PDFLIB 6.4.98
ccc         CALL PFTOPDG( BETA,DSQRT(Q2),XPQM)
         CALL PFTOPDG( BETA,Q2,XPQM)
         XPR=0.003D0
         ATMIN = (0.93827231D0*XPR)**2/(1.D0-XPR)
         DP = T_INT(ATMIN,1.D0,B0POM+2.D0*APPOM*LOG(1.D0/XPR)) *XPR**(-
     +   PP)
         DM = T_INT(ATMIN,1.D0,B0MES+2.D0*APMES*LOG(1.D0/XPR)) *XPR**(-
     +   PM)
         FLUX_P = 0.D0
         FLUX_M = 0.D0
         CALL FLUXH1(X_POM,T2,APPOM,B0POM,CP,PP,FLUX_P)
         CALL FLUXH1(X_POM,T2,APMES,B0MES,CM,PM,FLUX_M)
         FLUX_P = FLUX_P/(DP*XPR)
         FLUX_M = FLUX_M/(DM*XPR)
      ENDIF

      DO I=-6,6
         XPQ(I)=0.
         IF(ABS(I).LT.6) THEN
            IF(NPOM.EQ.-10 .OR. NPOM.EQ.-3) THEN
               XPQ(I)=REAL(XPQP(I)*FLUX_P)
            ELSEIF(NPOM.EQ.-11 ) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M)
            ELSEIF(NPOM.EQ.-12) THEN
               XPQ(I)=REAL(XPQM(I)*FLUX_M+XPQP(I)*FLUX_P)
            ENDIF
            IF(XPQ(I).LT.1E-10) XPQ(I)=0.
         ENDIF
         IF(FLAG) WRITE(6,*)'I,XPQ(I)',I,XPQ(I)
      ENDDO

      FIRST=.FALSE.

      RETURN
      END
