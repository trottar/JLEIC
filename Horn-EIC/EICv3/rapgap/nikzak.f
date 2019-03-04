*CMZ :  2.01/09 05/03/96  11.14.58  by  Hannes Jung
*CMZ :  2.01/04 24/01/96  16.10.26  by  Hannes Jung
*CMZ :  2.01/03 20/01/96  17.11.46  by  Hannes Jung
*CMZ :  2.01/01 05/01/96  12.53.32  by  Hannes Jung
*CMZ :  2.01/00 01/01/96  17.20.42  by  Hannes Jung
*CMZ :  2.00/30 31/12/95  15.10.06  by  Hannes Jung
*-- Author :    Hannes Jung   31/12/95

      SUBROUTINE NIKZAK(XBJ,XP,Q2,ABS_T,
     +           D2S_VL,D2S_VS,D2S_VC,D2S_SL,D2S_SS,D2S_SC)

*************************************************************************
*
* author: M.Grothe , 5.1.96
*
* D^2 SIGMA/(DT DXP) (GAMMA*P->XP) ACCORDING TO
* Genovese, Nikolaev, Zakharov, Cern/Th-95/13
*
* input : kinematic variables x_bj, x_pomeron, Q_sqrt, ABS(t)
* output: D^2 SIGMA/(DT DX_pomeron) in mbarn/GeV^2
*         D2S_VL  valence part with light quarks in final state
*         D2S_VS  valence part with strange quarks in final state
*         D2S_VC  valence part with charm quarks in final state
*         D2S_SL  sea part with light quarks in final state
*         D2S_SS  sea part with strange quarks in final state
*         D2S_SC  sea part with charm quarks in final state
*
*
*************************************************************************

      IMPLICIT NONE

      REAL XBJ,XP,Q2,ABS_T,BETA
      REAL CONST_VAL,CONST_SEA,D2S_VAL,D2S_SEA,DENOM_VAL,DENOM_SEA
      REAL D2S_VL,D2S_VS,D2S_VC,D2S_SL,D2S_SS,D2S_SC

      REAL ALPHA,PI

      REAL FLUX_SEA,FLUX_VAL,B_EL,B_3P,C_VAL,C_SEA,F_VAL,F_SEA
      REAL XP_0,SIGTOTPP

      PI        = 4.0*ATAN(1.0)
      ALPHA     = 1./137.

*!  in mbarn
      SIGTOTPP  = 40.
      XP_0      = 0.03

*!  in GeV^-2
      B_3P      = 6.
*!  in GeV^-2
      B_EL      = 12.

      C_VAL     = 0.27
      C_SEA     = 0.063

*******************************************************************************
C IN CASE XP < XP_0:
C  FLUX_VAL  = (XP_0/XP)**0.259 * (XP+0.00149)**0.2142 / (XP_0+0.00149)**0.2142
C            = CONST_VAL * (XP+0.00149)**0.2142 / XP**0.259
C  CONST_VAL = XP_0**0.259 / (XP_0+0.00149)**0.2142
C            = 0.845796
C  FLUX_SEA  = (XP_0/XP)**0.58 * (XP+0.0023)**0.48 / (XP_0+0.0023)**0.48
C            = CONST_SEA * (XP+0.0023)**0.48 / XP**0.58
C  CONST_SEA = XP_0**0.58 / (XP_0+0.0023)**0.48
C            = 0.679694
C ELSE:
C   FLUX_VAL  =1.
C   FLUX_SEA  =1.
C F_VAL     = C_VAL * BETA * (1-BETA)
C F_SEA     = C_SEA * (1-BETA)**2
*******************************************************************************

      CONST_VAL = 0.845796
      CONST_SEA = 0.679694

      IF (XP.GT.XP_0) THEN
         FLUX_VAL =1.
         FLUX_SEA =1
      ELSE
         FLUX_VAL = CONST_VAL * (XP+0.00149)**0.2142 / XP**0.259
         FLUX_SEA = CONST_SEA * (XP+0.0023)**0.48 / XP**0.58
      ENDIF

      BETA      = XBJ/XP

      F_VAL     = C_VAL * BETA * (1-BETA)
      F_SEA     = C_SEA * (1-BETA) * (1-BETA)

      D2S_VAL   = PI/4. * ALPHA * SIGTOTPP / Q2 * 1/XP
     +             * FLUX_VAL * F_VAL * EXP(-B_EL*ABS_T)

      D2S_SEA   = PI/4. * ALPHA * SIGTOTPP / Q2 * 1/XP
     +             * FLUX_SEA * F_SEA * EXP(-B_3P*ABS_T)
      IF(D2S_VAL.LE.0.0) D2S_VAL = 1.E-20
*******************************************************************************
C flavour decomposition for valence part:
C    A(u)=A(ubar)=A(d)=A(dbar)=0.20 of fraction in normal DIS
C    A(s)=A(sbar)=0.11
C    A(c)=A(cbar)=0.02
C flavour decomposition for sea part:
C    A(u)=A(ubar)=A(d)=A(dbar)=0.048 of fraction in normal DIS
C    A(s)=A(sbar)=0.040
C    A(c)=A(cbar)=0.009
*******************************************************************************

      DENOM_VAL  = 10./9.*0.20 + 2./9.*0.11 + 8./9.*0.02
      DENOM_SEA  = 10./9.*0.048 + 2./9.*0.040 + 8./9.*0.009

      D2S_VL     = D2S_VAL * 0.20*10./9./DENOM_VAL
      D2S_VS     = D2S_VAL * 0.11*2./9./DENOM_VAL
      D2S_VC     = D2S_VAL * 0.02*8./9./DENOM_VAL
      D2S_SL     = D2S_SEA * 0.048*10./9./DENOM_SEA
      D2S_SS     = D2S_SEA * 0.040*2./9./DENOM_SEA
      D2S_SC     = D2S_SEA * 0.009*8./9./DENOM_SEA
c      write(6,*) 'NIKZAK',XBJ,XP,Q2,ABS_T,
c     +           D2S_VL,D2S_VS,D2S_VC,D2S_SL,D2S_SS,D2S_SC
      RETURN
      END
