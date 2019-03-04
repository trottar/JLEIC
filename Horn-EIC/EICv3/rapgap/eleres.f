*CMZ :  2.08/04 22/12/99  15.39.23  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  15.59.10  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  19.15.01  by  Hannes Jung
*-- Author :    Hannes Jung   02/03/97
      SUBROUTINE ELERES(WT1)
C
C     Matrix elements for resolved photon processes
C     gg --> qq_bar
C     g + g --> g + g
C     g + q --> g + q
C     qq_bar --> gg
C     q + q_bar --> q + q_bar
C     qq --> qq
C
C
      IMPLICIT NONE
*KEEP,RGLUJETS.
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      REAL ULMASS
      DOUBLE PRECISION DOT1,DPLU,DLANGL
      EXTERNAL DLANGL,ULMASS,DPLU,LUCHGE,DOT1
C      SAVE

*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGPARA.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
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


*KEEP,RGLUCO.
      REAL PLEPIN,PPIN
      INTEGER KE,KP,KEB,KPH,KGL,KPA,NFRAG,ILEPTO,IFPS,IHF,IALMKT
      INTEGER INTER,ISEMIH
      INTEGER NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT,NFLAV,NFLQCDC
      COMMON/LUCO  /KE,KP,KEB,KPH,KGL,KPA,NFLAV,NFLQCDC
      COMMON/INPU  /PLEPIN,PPIN,NFRAG,ILEPTO,IFPS,IHF,IALMKT,INTER,
     +              ISEMIH
      COMMON/HARD/ NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT
      INTEGER IHFLA
      COMMON/HFLAV/ IHFLA
C      SAVE

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGLUDAT2.
      REAL PMAS,PARF,VCKM
      INTEGER KCHG
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C      SAVE

*KEEP,RGSTRU.
      REAL SCAL1,XPD1,SCAL2,XPD2
      COMMON/STRU/ SCAL1,XPD1,SCAL2,XPD2
C     SAVE
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
	Integer ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
	Integer IDEBUG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON/ INTERN/IDEBUG

	Integer IRFUL
      REAL XPQ(-6:6),XGQ(-6:6)
      REAL XGQ1,XPQ2,XGAM,XPRO,XGQ1_HF,XPQ2_HF
      REAL SNGL
c IRFUL=1 using full ME's
c IRFUL=0 using ME's from Bengtsson
      PARAMETER (IRFUL=0)
      DOUBLE PRECISION ZETA3
C---ZETA3=RIEMANN ZETA FUNCTION(3)
      PARAMETER (ZETA3=1.202056903159594D0)
C---paramtric fit to the function used in the Bartels/Forshaw/et al.
calculation
	Double Precision wt1,ALPHA_S,ALPHAS,SH,UH,TH,VTH,VUH
	Double Precision SUMA,SUMB,SUMC,SUMC1,SUMC1T,SUMD,SUMD1,SUMD2
	Double Precision SUMF,SUME,SUME1,SUME2,SUMTE
	Double Precision SUMTF,SUMF1,SUMF2,SUMG,SUMG1,SUTO,SUM2
	Double Precision A1,A2,A3,ASUM,FLSUM,FLRN
	Double Precision QTOT,QSUM,CF,CA,yh,omeg0,delta
	Double Precision RNTEST
	Integer I,NFG,NFLG,NFLP,NFLT,NFLTG,NFLTP,NFLTPC,KPF
	Double Precision draprn
      DOUBLE PRECISION alphasf,gmass
      PARAMETER (alphasf=0.23D0)
      PARAMETER (gmass=.8D0)
      wt1 = 0.d0
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      if(xpr.GT.1.d0) write(6,*) ' eleres xp2>1',xpr,ipro
      IF(IDIR.EQ.1) THEN
         XPRO = XPR
      ELSE
         XPRO = XPR/XFGKI
      ENDIF
      if(xpro.ge.0.999) return
      CALL PYSTFU(KINT(2,2),XPRO,SNGL(Q2Q),XPQ)
c call virtual photon structure function
      XGAM = XEL/YY
      if(xgam.gt.1.0) then
         write(6,*) ' eleres xpr,x_gam,yy',xpr,xgam,yy,q2q,q2
      endif
      if(xgam.gt.0.999) return
      CALL RYSTGA(XGAM,SNGL(Q2Q),SNGL(Q2),XGQ)
      if(xgq(0).ne.xgq(0)) then
         write(6,*) ' eleres x_gam,Q2Q,Q2,xgq(0)',xgam,q2q,q2,xgq(0)
      endif
c      write(6,*) ' eleres ',xgq
      XPQ2 = 0.
      XGQ1 = 0.
      DO 10 I=1,NFLAV
         XPQ2 = XPQ2+XPQ(I) + XPQ(-I)
         XGQ1 = XGQ1+XGQ(I) + XGQ(-I)
   10 CONTINUE
      IF(IHFLA.GE.4) THEN
         XPQ2_HF = XPQ(IHFLA) + XPQ(-IHFLA)
         XGQ1_HF = XGQ(IHFLA) + XGQ(-IHFLA)
      ENDIF
      SH = 2.D0 * AM(1)**2 + 2.D0 * DOT1(NF1,NF2)
      TH = AM(1)**2 - 2.D0 * DOT1(NIA1,NF1) + dot1(nia1,nia1)
      UH = AM(2)**2 - 2.D0 * DOT1(NIA1,NF2) + dot1(nia1,nia1)
      VTH = - 2.D0 * DOT1(NIA1,NF1)
      VUH = - 2.D0 * DOT1(NIA1,NF2)
      SUMA = 0.D0
      SUMB = 0.D0
      SUMC = 0.D0
      SUMD = 0.D0
      SUME = 0.D0
      SUMF = 0.D0
      SUMG = 0.D0
      IF(IRPA.NE.0) THEN
         IF(IHFLA.LT.4) THEN
C.....MATRIX ELEMENT for gg --> qq_barp(466)
            IF(IRFUL.EQ.1) THEN
               SUMA = 2.d0*th*uh/th**2/12.d0 + 2.d0*th*uh/uh**2/12.d0 -
     +         3.d0/16.d0*4.d0*(1.d0 -uh*th/sh**2) + 3.d0/32.d0*4.d0
            ELSE
c this was from bengtsson
               SUMA = uh/th -2.d0*uh**2/sh**2 + th/uh - 2.d0*th**2/sh**
     +         2
               SUMA = SUMA/6.d0
            ENDIF
c include sum over NFG flavors and parton densities
            NFG = 3
c         IF(SH.GT.4.D0*DBLE(PMAS(4,1)**2)) NFG = 4
c         IF(SH.GT.4.D0*DBLE(PMAS(5,1)**2)) NFG = 5

            SUMA = SUMA * DFLOAT(NFG) * DBLE(XGQ(0)) * DBLE(XPQ(0))
         ELSEIF(IHFLA.GE.4) THEN
c heavy flavour gg --> QQ_bar
            SUMA = (32.D0/vth/vuh - 8d0*9D0/sh**2)*
     +   (vth**2/4d0 + vuh**2/4d0 + AM(1)**2*sh - AM(1)**2*sh/vth/vuh)
            SUMA = SUMA/16d0/3d0
            SUMA = SUMA  * DBLE(XGQ(0)) * DBLE(XPQ(0))
         ENDIF
      ENDIF
      IF(IRPB.NE.0) THEN
         IF(IRFUL.EQ.1) THEN
C.....MATRIX ELEMENT for gg --> gg
            SUMB = 9.d0/8.d0*(17.d0/2.d0 - 4.d0*uh*sh/th**2 + 17.d0/
     +      2.d0 - 4.d0*sh*th/uh**2 + 17.d0/2.d0 - 4.d0*uh*th/sh**2 +
     +      27.d0) + 9.d0/16.d0*(15.d0 - sh**2/th/uh + 15.d0 - uh**2/
     +      th/sh - 15.d0 + th**2/sh/uh) - 3.d0*9.d0/8.d0*81.d0/4.d0
         ELSE
c this is from bengtsson
            SUMB = 9.d0/4.d0*( sh**2/th**2+2.d0*sh/th+3.d0+2.d0*th/sh+
     +      th** 2/sh**2 +uh**2/sh**2+2.d0*uh/sh+3.d0+2.d0*sh/uh+sh**2/
     +      uh**2 +th**2/uh**2+2.d0*th/uh+3.d0+2.d0*uh/th+uh**2/th**2 )
     +      * 0.5d0
ccc         write(6,*) sumb,xgq(0),xpq(0),uh,th,sh
         ENDIF
         SUMB = SUMB * DBLE(XGQ(0)) * DBLE(XPQ(0))

      ENDIF
      IF(IRPC.NE.0) THEN
C.....MATRIX ELEMENT for qg --> qg
         IF(IRFUL.EQ.1) THEN
            SUMC1 = 2.d0*(1.d0-uh*sh/th**2)-4.d0/9.d0*(sh/uh+uh/sh)-
     +      1.d0
c  this was from bengtsson
         ELSE
            SUMC1 = 4.d0/9.d0*(2.d0*uh**2/th**2-uh/sh+2.d0*sh**2/th**2-
     +      sh/ uh)
         ENDIF
         sumc1t = SUMC1 *(DBLE(XGQ1) * DBLE(XPQ(0)))
c      write(6,*) xgq1,xgg1,xpq2,xpg2
         SUMC = SUMC1 *(DBLE(XGQ1)*DBLE(XPQ(0))+DBLE(XGQ(0))*DBLE(XPQ2)
     +   )
         IF(IHFLA.GE.4) THEN
            SUMC=0.D0
            IF(sh.ge.4.*ULMASS(IHFLA)) THEN
               SUMC = SUMC1* (DBLE(XGQ1_HF)*DBLE(XPQ(0))+DBLE(XGQ(0))*
     +         DBLE(XPQ2_HF))
            ENDIF
         ENDIF
      ENDIF
      IF(IRPD.NE.0) THEN
C.....MATRIX ELEMENT for qq_bar --> gg
         IF(IRFUL.EQ.1) THEN
            SUMD1 = 16.d0/27.d0*2.d0*uh/th + 16.d0/27.d0*2d0*th/uh
     +      -4.d0/3.d0*4.d0*(1.d0-uh*th/sh**2) +2.d0/3.d0*4.d0
            SUMD2 = 0.D0
         ELSE
            SUMD1 =32.d0/27.d0*(uh/th-2.d0*uh**2/sh**2)
            SUMD2 =32.d0/27.d0*(th/uh-2.d0*th**2/sh**2)
         ENDIF
c
         SUMD = (SUMD1+SUMD2) * DBLE(XGQ(1)*XPQ(-1)+ XGQ(2)*XPQ(-2)+
     +   XGQ(3)*XPQ(-3)+ XGQ(4)*XPQ(-4)+ XGQ(-1)*XPQ(1)+ XGQ(-2)*XPQ(2)
     +   + XGQ(-3)*XPQ(3)+ XGQ(-4)*XPQ(4))
      ENDIF
      IF(IRPE.NE.0) THEN
C.....MATRIX ELEMENT for qq_bar --> qq_bar

         IF(IHFLA.LT.4) THEN

            SUME1 = 4.d0/9.d0*(sh**2 + uh**2)/th**2
            SUME2 = 4.d0/9.d0*(th**2 + uh**2)/sh**2
         ELSEIF(IHFLA.GE.4) THEN
c heavy flavour qq --> QQ_bar
            SUME1 = 0.D0
            SUME2=4.d0/9.d0*
     +      ((vth**2 + vuh**2)/sh**2 + 4d0*AM(1)**2/sh**2)
         ENDIF
c
c first take t channel: q_i + q_bar_j --> q_i + q_bar_j
         SUME1 = SUME1 * (DBLE(XGQ(1))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+
     +   XPQ(-1)) + DBLE(XGQ(2))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1))
     +   + DBLE(XGQ(3))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1)) +
     +   DBLE(XGQ(4))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1)) +
     +   DBLE(XGQ(-1))*
     +   DBLE(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(-2))*DBLE(XPQ(4)+
     +   XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(-3))*DBLE(XPQ(4)+XPQ(3)+
     +   XPQ(2)+XPQ(1)) + DBLE(XGQ(-4))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+
     +   XPQ(1)))
c now s channel q_i + q_bar_i --> q_k + q_bar_k
         SUME2 = SUME2 * DBLE(XGQ(1)*XPQ(-1)+ XGQ(2)*XPQ(-2)+ XGQ(3)*
     +   XPQ(-3)+ XGQ(4)*XPQ(-4)+ XGQ(-1)*XPQ(1)+ XGQ(-2)*XPQ(2)+
     +   XGQ(-3)*XPQ(3)+ XGQ(-4)*XPQ(4))
         SUME = SUME1 + SUME2
      ENDIF
      IF(IRPF.NE.0) THEN
C.....MATRIX ELEMENT for qq --> qq, q_bar q_bar --> q_bar q_bar


         SUMF1 = 4.d0/9.d0*(sh**2 + uh**2)/th**2
         SUMF2 = 4.d0/9.d0*(sh**2 + th**2)/uh**2
c

         SUMF1 = SUMF1 * (DBLE(XGQ(1))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)
     +   ) + DBLE(XGQ(2))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(
     +   3))*DBLE(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(4))*
     +   DBLE(XPQ(4)+XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(-1))*DBLE(XPQ(-4)
     +   +XPQ(-3)+XPQ(-2)+XPQ(-1)) + DBLE(XGQ(-2))*DBLE(XPQ(-4)+XPQ(-3)
     +   +XPQ(-2)+XPQ(-1)) + DBLE(XGQ(-3))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)
     +   +XPQ(-1)) + DBLE(XGQ(-4))*DBLE(XPQ(-4)+XPQ(-3)+XPQ(-2)+XPQ(-1)
     +   ))
         SUMF2 = SUMF2 * (DBLE(XGQ(1))*DBLE(XPQ(1)) + DBLE(XGQ(2))*
     +   DBLE(XPQ(2)) + DBLE(XGQ(3))*DBLE(XPQ(3)) + DBLE(XGQ(4))*
     +   DBLE(XPQ(4)) + DBLE(XGQ(-1))*DBLE(XPQ(-1)) + DBLE(XGQ(-2))*
     +   DBLE(XPQ(-2)) + DBLE(XGQ(-3))*DBLE(XPQ(-3)) + DBLE(XGQ(-4))*
     +   DBLE(XPQ(-4)))
         SUMF = SUMF1 + SUMF2
      ENDIF
      IF(IRPG.NE.0) THEN
c cross sections from B. Cox J. Forshaw hep-ph/9805206
         SUMG1 = 0.D0
         if(IRPG.eq.1) then
c eq. 3 Mueller tang
            CF = 4.D0/3.D0
            CA = 3.D0
            yh=LOG(SH/(-TH))+1.D0
            omeg0 = CA*alphasf/pi*4*dlog(2.d0)
            SUMG1=16.D0*PI*2*PI**3*(CF*alphasf)**4
     +    *(SH/TH)**2 *dexp(2.d0*omeg0*yh)/(7.D0*alphasf*CA*zeta3*yh)**3
         elseif(IRPG.eq.2) then
c eq 4. massive gluon
            delta=sqrt(1-((4*gmass**2)/TH))
            SUMG1=(SH/TH)**2/delta**2*(log((delta+1)/(delta-1)))**2*
     +      PI**2*(alphasf)**4
         endif
         SUMG=SUMG1*DBLE(XPQ2)*DBLE(XGQ1)
      ENDIF
c now the total
      SUTO = SUMA + SUMB + SUMC + SUMD + SUME + SUMF + SUMG
c      SUTO = SUMA + SUMB + SUMC + SUMD + SUME + SUMF
      SUM2 = SUTO*16.D0*PI*PI*ALPHA_S*ALPHA_S
      WT1=SUM2
      RNTEST = draprn()
      IF((SUMA+SUMB+SUMC+SUMD+SUME+SUMF)/SUTO.LT.RNTEST) THEN
         IRESPRO = 7
      ELSEIF((SUMA+SUMB+SUMC+SUMD+SUME)/SUTO.LT.RNTEST) THEN
         IRESPRO = 6
      ELSEIF((SUMA+SUMB+SUMC+SUMD)/SUTO.LT.RNTEST) THEN
         IRESPRO = 5
      ELSEIF((SUMA+SUMB+SUMC)/SUTO.LT.RNTEST) THEN
         IRESPRO = 4
      ELSEIF((SUMA+SUMB)/SUTO.LT.RNTEST) THEN
         IRESPRO = 3
      ELSEIF(SUMA/SUTO.LT.RNTEST) THEN
         IRESPRO = 2
      ELSE
         IRESPRO = 1
      ENDIF
      ICOLORA = 0
c      write(6,*) ' eleres ',IRESPRO
      IF(IRESPRO.EQ.1) THEN
c get color configuartion for gg --> qq_bar
         A1 = uh/th - 2.d0*uh**2/sh**2
         A2 = th/uh - 2.d0*th**2/uh**2
         ASUM = A1 + A2
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
         IF(IHFLA.LT.4) THEN
c select flavour of outgoing q q_bar
            FLRN = draprn()
            FLSUM = 1.D0/3.D0
            IF(FLRN.LT.FLSUM) THEN
               KPF = 1
            ELSEIF(FLRN.LT.2.D0*FLSUM) THEN
               KPF = 2
            ELSEIF(FLRN.LE.3.D0*FLSUM) THEN
               KPF = 3
            ENDIF
         ELSE
            KPF = IHFLA
         ENDIF
         IF(draprn().LE.0.5) THEN
            K(NF1,2) = KPF
            K(NF2,2) = - KPF
         ELSE
            K(NF1,2) = - KPF
            K(NF2,2) = KPF
         ENDIF

      ELSEIF(IRESPRO.EQ.2) THEN
c get color configuartion for gg --> gg
         K(NF1,2)=21
         K(NF2,2)=21
         A1 = sh**2/th**2+2.d0*sh/th+3.d0+2.d0*th/sh+th**2/sh**2
         A2 = uh**2/sh**2+2.d0*uh/sh+3.d0+2.d0*sh/uh+sh**2/uh**2
         A3 = th**2/uh**2+2.d0*th/uh+3.d0+2.d0*uh/th+uh**2/th**2
         ASUM = A1+A2+A3
         RNTEST = draprn()
         ICOLORA = 0
         IF((A1+A2)/ASUM.lt.rntest) THEN
            ICOLORA=3
         ELSEIF(A1/ASUM.lt.rntest) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
      ELSEIF(IRESPRO.EQ.3) THEN
c where is the quark ( from gamma or proton)
         IF(SUMC1t/SUMC.GT.draprn()) THEN
            IF(IHFLA.LT.4) THEN
               NFLT = -NFLAV-1
               QSUM = -DBLE(XGQ1)*draprn()
   20          NFLT = NFLT + 1
               IF(NFLT.EQ.0) GOTO 20
c         write(6,*) ' q of gamma ',QSUM,XGQ1,NFLT,XGQ(NFLT)
               QSUM = QSUM + DBLE(XGQ(NFLT))
               IF(QSUM.LE.0.0D0) GOTO 20
               IF(NFLT.GT.NFLAV) write(6,*) ' eleres NFL > NFLAV ',
     +         NFLT, NFLAV
            ELSE
               NFLT = IHFLA
            ENDIF
            IF(draprn().GT.0.5) NFLT=-IHFLA
            K(NIA1,2) = NFLT
            K(NIR1,2) = -NFLT
            K(NIA2,2) = 21
            K(NF1,2) = NFLT
            K(NF2,2) = 21
         ELSE
            IF(IHFLA.LT.4) THEN
               NFLT = -NFLAV-1
               QSUM = -DBLE(XPQ2)*draprn()
   30          NFLT = NFLT + 1
               IF(NFLT.EQ.0) GOTO 30
c         write(6,*) ' q of proton ',QSUM,XPQ2,NFLT,XPQ(NFLT)
               QSUM = QSUM + DBLE(XPQ(NFLT))
               IF(QSUM.LE.0.0D0) GOTO 30
               IF(NFLT.GT.NFLAV) write(6,*) ' eleres NFL > NFLAV ',
     +         NFLT, NFLAV
            ELSE
               NFLT = IHFLA
            ENDIF
            IF(draprn().GT.0.5) NFLT=-IHFLA
            K(NIA1,2) = 21
            K(NIR1,2) = 21
            K(NIA2,2) = NFLT
            K(NF1,2) = 21
            K(NF2,2) = NFLT
         ENDIF
c         call dulist(1)
c get color configuartion for gq --> gq
         A1 = 2.d0*uh**2/th**2 - uh/sh
         A2 = 2.d0*sh**2/th**2 - sh/th
         ASUM = A1 + A2
         ICOLORA = 0
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
      ELSEIF(IRESPRO.EQ.4) THEN
c get color configuartion for gg --> qq_bar
         A1 = uh/th - 2.d0*uh**2/sh**2
         A2 = th/uh - 2.d0*th**2/uh**2
         ASUM = A1 + A2
         IF(A1/ASUM.lt.draprn()) THEN
            ICOLORA=2
         ELSE
            ICOLORA=1
         ENDIF
         IF(IHFLA.LT.4) THEN
            QTOT=DBLE(XGQ(1)*XPQ(-1)+ XGQ(2)*XPQ(-2)+ XGQ(3)*XPQ(-3)+
     +      XGQ(-1)*XPQ(1)+ XGQ(-2)*XPQ(2)+ XGQ(-3)*XPQ(3))
            NFLT = -3-1
            QSUM = -QTOT*draprn()
   40       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 40
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 40
            IF(NFLT.GT.NFLAV) write(6,*) ' eleres NFL > NFLAV ',NFLT,
     +      NFLAV
         ELSE
            NFLT = IHFLA
            IF(draprn().GT.0.5) NFLT= -IHFLA
         ENDIF
c         write(6,*) ' IRPE 2nd ',NFLT
         K(NIA1,2) = -NFLT
         K(NIR1,2) = NFLT
         K(NIA2,2) = NFLT
         K(NF1,2) = 21
         K(NF2,2) = 21
      ELSEIF(IRESPRO.EQ.5) THEN
c qq_bar --> qq_bar
         SUMTE = SUME1+SUME2
         IF(SUME1/SUMTE.lt.draprn()) THEN
            ICOLORA = 2
            QTOT  =     DBLE(XGQ(1)*XPQ(-1)+
     +                    XGQ(2)*XPQ(-2)+
     +                    XGQ(3)*XPQ(-3)+
     +                    XGQ(-1)*XPQ(1)+
     +                    XGQ(-2)*XPQ(2)+
     +                    XGQ(-3)*XPQ(3))
            NFLT = -3-1
            QSUM = -QTOT*draprn()
   50       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 50
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(-NFLT))
            IF(QSUM.LE.0.0D0) GOTO 50
            IF(NFLT.GT.NFLAV) write(6,*) ' eleres NFL > NFLAV ',NFLT,
     +      NFLAV
            K(NIA1,2) = -NFLT
            K(NIR1,2) = NFLT
            K(NIA2,2) = NFLT
c select flavour of outgoing q q_bar
            IF(IHFLA.LT.4) THEN
               FLRN = draprn()
               FLSUM = 1.D0/3.D0
               IF(FLRN.LT.FLSUM) THEN
                  KPF = 1
               ELSEIF(FLRN.LT.2.D0*FLSUM) THEN
                  KPF = 2
               ELSEIF(FLRN.LE.3.D0*FLSUM) THEN
                  KPF = 3
               ENDIF
            ELSE
               KPF = IHFLA
            ENDIF
            IF(K(NIA1,2).LE.0) THEN
               K(NF1,2) = -KPF
               K(NF2,2) =  KPF
            ELSE
               K(NF1,2) = KPF
               K(NF2,2) = -KPF
            ENDIF
         ELSE
            ICOLORA = 1
            QTOT = (DBLE(XGQ(1))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)) +
     +                DBLE(XGQ(2))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)) +
     +                DBLE(XGQ(3))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)) +
     +                DBLE(XGQ(-1))*DBLE(XPQ(3)+XPQ(2)+XPQ(1)) +
     +                DBLE(XGQ(-2))*DBLE(XPQ(3)+XPQ(2)+XPQ(1)) +
     +                DBLE(XGQ(-3))*DBLE(XPQ(3)+XPQ(2)+XPQ(1)))
            NFLTG = -3-1
            QSUM = -QTOT*draprn()
   60       NFLTG = NFLTG + 1
            IF(NFLTG.EQ.0) GOTO 60
            IF(NFLTG.LE.0) THEN
               NFLTP = 0
               NFLTPC = 3
            ELSE
               NFLTP = -3-1
               NFLTPC = 0
            ENDIF
   70       NFLTP = NFLTP + 1
            IF(NFLTP.EQ.0) GOTO 70
            IF(NFLTP.GT.NFLTPC) GOTO 60
            QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
            IF(QSUM.LE.0.D0.AND.NFLTP.LE.NFLTPC) GOTO 70
            IF(QSUM.LE.0.0D0) GOTO 60
            IF(NFLTG.GT.NFLAV) write(6,*) ' eleres NFL_gam > NFLAV ',
     +      NFLTG, NFLAV
            IF(NFLTP.GT.NFLAV) write(6,*) ' eleres NFL_gp > NFLAV ',
     +      NFLTP, NFLAV
            K(NIA1,2) = NFLTG
            K(NIR1,2) = -NFLTG
            K(NIA2,2) = NFLTP
            K(NF1,2) = NFLTG
            K(NF2,2) = NFLTP
         ENDIF
      ELSEIF(IRESPRO.EQ.6) THEN
C qq --> qq, q_bar q_bar --> q_bar q_bar
         SUMTF = SUMF1 + SUMF2
         IF(SUMF1/SUMTF.lt.draprn()) THEN
            ICOLORA = 2
            QTOT = (DBLE(XGQ(1))*DBLE(XPQ(1)) + DBLE(XGQ(2))*
     +      DBLE(XPQ(2)) + DBLE(XGQ(3))*DBLE(XPQ(3)) + DBLE(XGQ(-1))*
     +      DBLE(XPQ(-1)) + DBLE(XGQ(-2))*DBLE(XPQ(-2)) + DBLE(XGQ(-3))
     +      *DBLE(XPQ(-3)))
            NFLT = -3-1
            QSUM = -QTOT*draprn()
   80       NFLT = NFLT + 1
            IF(NFLT.EQ.0) GOTO 80
            QSUM = QSUM + DBLE(XPQ(NFLT)*XGQ(NFLT))
            IF(QSUM.LE.0.0D0) GOTO 80
            IF(NFLT.GT.NFLAV) write(6,*) ' eleres NFL > NFLAV ',NFLT,
     +      NFLAV
            K(NIA1,2) = NFLT
            K(NIR1,2) = -NFLT
            K(NIA2,2) = NFLT
            K(NF1,2) = NFLT
            K(NF2,2) = NFLT
         ELSE
            ICOLORA = 1
            QTOT = (DBLE(XGQ(1))*DBLE(XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(
     +      2))*DBLE(XPQ(3)+XPQ(2)+XPQ(1)) + DBLE(XGQ(3))*DBLE(XPQ(3)+
     +      XPQ(2)+XPQ(1)) + DBLE(XGQ(-1))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)
     +      ) + DBLE(XGQ(-2))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)) + DBLE(XGQ(
     +      -3))*DBLE(XPQ(-3)+XPQ(-2)+XPQ(-1)))
            NFLTG = -3-1
            QSUM = -QTOT*draprn()
   90       NFLTG = NFLTG + 1
            IF(NFLTG.EQ.0) GOTO 90
            IF(NFLTG.LE.0) THEN
               NFLTP = -3-1
               NFLTPC = 0
            ELSE
               NFLTP = 0
               NFLTPC = 3
            ENDIF
  100       NFLTP = NFLTP + 1
            IF(NFLTP.GT.NFLTPC) GOTO 90
            IF(NFLTP.EQ.0) GOTO 100
            QSUM = QSUM + DBLE(XPQ(NFLTP)*XGQ(NFLTG))
            IF(QSUM.LE.0.D0.AND.NFLTP.LE.NFLTPC) GOTO 100
            IF(QSUM.LE.0.0D0) GOTO 90
            IF(NFLTG.GT.NFLAV) write(6,*) ' eleres NFL_gam > NFLAV ',
     +      NFLTG, NFLAV
            IF(NFLTP.GT.NFLAV) write(6,*) ' eleres NFL_gp > NFLAV ',
     +      NFLTP, NFLAV
            K(NIA1,2) = NFLTG
            K(NIR1,2) = -NFLTG
            K(NIA2,2) = NFLTP
            K(NF1,2) = NFLTG
            K(NF2,2) = NFLTP
         ENDIF
      ELSEIF(IRESPRO.EQ.7) THEN
C qq --> qq, q_bar q_bar --> q_bar q_bar color singlet exchange
         NFLP = -NFLAV-1
         QSUM = -DBLE(XPQ2)*draprn()
  110    NFLP = NFLP + 1
         IF(NFLP.EQ.0) GOTO 110
c         write(6,*) ' q of proton ',QSUM,XPQ2,NFLP,XPQ(NFLP)
         QSUM = QSUM + DBLE(XPQ(NFLP))
         IF(QSUM.LE.0.0D0) GOTO 110
         IF(NFLP.GT.NFLAV) write(6,*) ' eleres NFLP > NFLAV ',NFLP,
     +   NFLAV
         NFLG = -NFLAV-1
         QSUM = -DBLE(XGQ1)*draprn()
  120    NFLG = NFLG + 1
         IF(NFLG.EQ.0) GOTO 120
c         write(6,*) ' q of gamma ',QSUM,XGQ1,NFLG,XGQ(NFLG)
         QSUM = QSUM + DBLE(XGQ(NFLG))
         IF(QSUM.LE.0.0D0) GOTO 120
         IF(NFLG.GT.NFLAV) write(6,*) ' eleres NFLG > NFLAV ',NFLG,
     +   NFLAV

         K(NIA1,2) = NFLG
         K(NIR1,2) = -NFLG
         K(NIA2,2) = NFLP
         K(NF1,2) = NFLG
         K(NF2,2) = NFLP
      ENDIF
      SCAL1 = SNGL(Q2Q)
      XPD1 = XGQ(K(NIA1,2))
      IF(WT1.LE.0.D0.AND.IDEBUG.EQ.1) THEN
         write(6,*) ' eleres wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         write(6,*) ' eleres xgq = ',(xgq(i),i=-6,6)
         write(6,*) ' eleres xpq = ',(xpq(i),i=-6,6)
      ENDIF
      IF(WT1.LT.0.D0) THEN
         write(6,*) ' eleres wt1 = ',wt1,' Q2 = ',Q2,' Q2Q = ',Q2Q
         WT1 = 0.D0
      ENDIF
      RETURN
      END
