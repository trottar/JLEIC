*CMZ :  2.06/34 20/05/98  09.12.56  by  Hannes Jung
*CMZ :  2.06/30 18/03/98  13.30.04  by  Hannes Jung
*CMZ :  2.06/29 15/03/98  14.20.08  by  Hannes Jung
*CMZ :  2.06/19 19/11/97  10.41.40  by  Hannes Jung
*CMZ :  2.06/16 29/10/97  09.37.37  by  Hannes Jung
*CMZ :  2.06/15 22/10/97  17.43.35  by  Hannes Jung
*CMZ :  2.06/11 17/10/97  13.05.17  by  Hannes Jung
*CMZ :  2.06/02 22/08/97  03.04.20  by  Hannes Jung
*CMZ :  2.06/00 16/07/97  15.21.49  by  Hannes Jung
*CMZ :  2.05/21 26/06/97  17.47.14  by  Hannes Jung
*CMZ :  2.05/13 26/05/97  12.30.13  by  Hannes Jung
*CMZ :  2.05/09 06/05/97  18.38.07  by  Hannes Jung
*CMZ :  2.05/07 24/04/97  17.59.59  by  Hannes Jung
*CMZ :  2.05/06 23/04/97  06.59.50  by  Hannes Jung
*CMZ :  2.05/05 20/03/97  22.41.09  by  Hannes Jung
*CMZ :  2.05/00 05/03/97  08.10.05  by  Hannes Jung
*-- Author :    Hannes Jung   03/03/97
      SUBROUTINE RYSTGA(XS,SCAL,Q2S,XPG)
      REAL XS,SCAL,Q2S,XPG(-6:6)
      DOUBLE PRECISION X,Q2,P2,UGAM,DGAM,SGAM,GGAM
      DOUBLE PRECISION PP2
      DIMENSION XPDFGM(-6:6)
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEEP,RGPARAS.
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE

*KEEP,RGSCQ2.
      DOUBLE PRECISION SCALQ2
	Integer ISET,IP2
      COMMON/RESGAM/ISET,IP2,SCALQ2
*KEND.
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      CHARACTER*20 PARM(20)
      LOGICAL FIRST,PDFFIRST,TEST
      COMMON /W50516/PDFFIRST
      DATA NPDF/0/,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU/8*0.D0/
      DATA FIRST/.TRUE./
      DATA TEST/.FALSE./
c      DATA IGAMPDF/1/
      IGAMPDF=MSTP(56)
c IGAMPDF = 1  --> GRS - LO - VIRTUAL PHOTON PARAMETRIZATIONS
c IGAMPDF = 2  --> SaSgam version 2 - parton distributions of the photon
      X=DBLE(XS)
      Q2=DBLE(SCAL)
      P2=DBLE(Q2S)
      DO 10 I=-6,6
   10 XPG(I)=1.e-21
ccccc  10 XPG(I)=0.0

      IF(FIRST) THEN
         FIRST=.FALSE.
         write(6,*) ' parton densities in the photon '
         write(6,10000) SCALQ2
10000    FORMAT('  cut scale > ',F6.2,' Q2 ')
         IF(IGAMPDF.EQ.1) THEN
            write(6,*) ' GRS - LO - parametrisation'
         ELSEIF(IGAMPDF.EQ.2) THEN
            write(6,*) ' SaSgam version 2 (Schuler - Sjostrand) '
            write(6,*) ' ISET = ',ISET,' IPS = ',IP2
         ELSEIF(IGAMPDF.GT.1000) THEN
            write(6,*) ' PDFLIB used with ',IGAMPDF
            IF(OMEG2.GT.0D0) THEN
               write(6,*) ' virtual photon suppression a la Drees-'
     +         //'Godbole'
               write(6,*) ' L=(ln(scal+omega^2)/(Q^2+omega^2))/',
     +         '(ln(scal+omega^2)/omega^2)) '
               write(6,*) ' with omega**2 = ',OMEG2
            ENDIF
         ENDIF
      ENDIF
      IF(XS.GE.1.0) RETURN
c check that scale > gamma virtuality
      IF(Q2.LT.SCALQ2*1.0001*P2) RETURN
      IF(IGAMPDF.EQ.1) THEN
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*           G R S - LO - VIRTUAL PHOTON PARAMETRIZATIONS          *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                M. GLUECK, E.REYA, M. STRATMANN :                *
*                    PHYS. REV. D51 (1995) 3220                   *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.6   AND  5.E4                 *
*                       AND (!)  Q**2 > 5 P**2                    *
*        P**2 / GEV**2  BETWEEN   0.0   AND  10.                  *
*                       P**2 = 0  <=> REAL PHOTON                 *
*             X         BETWEEN  1.E-4  AND   1.                  *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*                                                                 *
*   PLEASE REPORT ANY STRANGE BEHAVIOUR TO :                      *
*          STRAT@HAL1.PHYSIK.UNI-DORTMUND.DE                      *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*    P2  = VIRTUALITY OF THE PHOTON IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION DIVIDED BY ALPHA_EM) :
*
********************************************************
         IF(Q2.GE.5D4) RETURN
         IF(Q2.LE.5.D0*P2.OR.P2.GT.10.D0) THEN

         ELSE
            call grspar(x,q2,p2,ugam,dgam,sgam,ggam)
            XPG(0)=SNGL(GGAM)/137.
            XPG(1)=SNGL(DGAM)/137.
            XPG(2)=SNGL(UGAM)/137.
            XPG(3)=SNGL(SGAM)/137.
            XPG(-1)=XPG(1)
            XPG(-2)=XPG(2)
            XPG(-3)=XPG(3)
            XPG(-4)=XPG(4)
            XPG(-5)=XPG(5)
            XPG(-6)=XPG(6)
            IF(TEST) THEN
               write(6,*) ' old GRS ',XPG(0),XPG(1)
               XXX=DBLE(XS)
               QQ=DBLE((MAX(0.,SCAL)))
C if you use PDFLIB version 3.  then uncoment the following 2 lines
C and comment the lines with 'cold pdf'
cold pdf  PARM(1)='MODE'
cold pdf  VAL(1)=MSTP(51)
C if you use PDFLIB 4.   then use the next lines.
CNEW
               ITEST = 3005004
               PARM(1) = 'NPTYPE'
               VAL(1) = INT(ITEST/1000000)
               PARM(2) = 'NGROUP'
               VAL(2) = DBLE(MOD(ITEST,1000000)/1000)
               PARM(3) = 'NSET'
               VAL(3) = DBLE(MOD(MOD(ITEST,1000000),1000))
CNEW
               NPDF=NPDF+1
               PDFFIRST=.FALSE.
               IF(NPDF.LE.1) PDFFIRST=.TRUE.
c call PDFSET each time, because when DIF with pi structure function is
c also selected one would get in confusion....
               CALL PDFSET(PARM,VAL)
cc         CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
c new calling sequence needed
               CALL STRUCTP(XXX,QQ,P2,IP2, UPV,DNV,USEA,DSEA,STR,CHM,
     +         BOT,TOP,GLU)
               Q2SUP = 1.D0
               XPG(0)=SNGL(GLU)*Q2SUP**2
               XPG(1)=SNGL(DSEA)*Q2SUP
               XPG(-1)=SNGL(DSEA)*Q2SUP
               XPG(2)=SNGL(USEA)*Q2SUP
               XPG(-2)=SNGL(USEA)*Q2SUP
               XPG(3)=SNGL(STR)*Q2SUP
               XPG(-3)=SNGL(STR)*Q2SUP
               XPG(4)=SNGL(CHM)*Q2SUP
               XPG(-4)=SNGL(CHM)*Q2SUP
               XPG(5)=SNGL(BOT)*Q2SUP
               XPG(-5)=SNGL(BOT)*Q2SUP
               XPG(6)=SNGL(TOP)*Q2SUP
               XPG(-6)=SNGL(TOP)*Q2SUP
               write(6,*) ' NEW GRS ',XPG(0),XPG(1)
            ENDIF
         ENDIF
      ELSEIF(IGAMPDF.EQ.2) THEN
         CALL SASGAM(ISET,XS,SCAL,Q2S,IP2,F2GM,XPDFGM)
C...The calling sequence is the following:
C     CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
C           IP2 : scheme used to evaluate off-shell anomalous component.
C               = 0 : recommended default, see = 7.
C               = 1 : dipole dampening by integration; very time-consuming.
C               = 2 : P_0^2 = max( Q_0^2, P^2 )
C               = 3 : P'_0^2 = Q_0^2 + P^2.
C               = 4 : P_{eff} that preserves momentum sum.
C               = 5 : P_{int} that preserves momentum and average
C                     evolution range.
C               = 6 : P_{eff}, matched to P_0 in P2 -> Q2 limit.
C               = 7 : P_{eff}, matched to P_0 in P2 -> Q2 limit.
C...The breakdown by component is stored in the commonblock SASCOM,
C               with elements as above.
C           XPVMD : rho, omega, phi VMD part only of output.
C           XPANL : d, u, s anomalous part only of output.
C           XPANH : c, b anomalous part only of output.
C           XPBEH : c, b Bethe-Heitler part only of output.
C           XPDIR : Cgamma (direct contribution) part only of output.
         DO 20 I=-6,6
   20    XPG(I)=XPDFGM(I)

      ELSEIF(IGAMPDF.GT.1000) THEN
C...Call PDFLIB structure functions.
         XXX=DBLE(XS)
         QQ=DBLE((MAX(0.,SCAL)))
         PARM(1) = 'NPTYPE'
         VAL(1) = INT(IGAMPDF/1000000)
         PARM(2) = 'NGROUP'
         VAL(2) = DBLE(MOD(IGAMPDF,1000000)/1000)
         PARM(3) = 'NSET'
         VAL(3) = DBLE(MOD(MOD(IGAMPDF,1000000),1000))
CNEW
c         call PDFSET each time, because when DIF with pi structure function is
c also selected one would get in confusion....
         NPDF=NPDF+1
         PDFFIRST=.FALSE.
         IF(NPDF.LE.1) PDFFIRST=.TRUE.
         CALL PDFSET(PARM,VAL)
c new calling sequence needed
c now here QQ is in GeV**2 and not as before in GeV
         IF(OMEG2.LT.0.D0) THEN
            Q2SUP = 1.D0
            PP2 = P2
         ELSE
            Q2SUP = SNGL(LOG((DBLE(SCAL)+OMEG2)/(P2+OMEG2))/ LOG((DBLE(
     +      SCAL)+OMEG2)/OMEG2))
            PP2 = 0.D0
         ENDIF
         CALL STRUCTP(XXX,QQ,PP2,IP2,
     +     UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
         XPG(0)=SNGL(GLU)*Q2SUP**2
         XPG(1)=SNGL(DSEA)*Q2SUP
         XPG(-1)=SNGL(DSEA)*Q2SUP
         XPG(2)=SNGL(USEA)*Q2SUP
         XPG(-2)=SNGL(USEA)*Q2SUP
         XPG(3)=SNGL(STR)*Q2SUP
         XPG(-3)=SNGL(STR)*Q2SUP
         XPG(4)=SNGL(CHM)*Q2SUP
         XPG(-4)=SNGL(CHM)*Q2SUP
         XPG(5)=SNGL(BOT)*Q2SUP
         XPG(-5)=SNGL(BOT)*Q2SUP
         XPG(6)=SNGL(TOP)*Q2SUP
         XPG(-6)=SNGL(TOP)*Q2SUP
      ELSE
         write(6,*) ' photon parametrisation not implemented ',IGAMPDF
      ENDIF
      RETURN
      END
