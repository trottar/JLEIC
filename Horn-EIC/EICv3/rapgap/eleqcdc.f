*CMZ :  2.08/00 06/06/99  15.49.36  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  18.23.37  by  Hannes Jung
*-- Author :    Hannes Jung   18/04/95
      SUBROUTINE ELEQCDC(WT1)
C    QCD Compton
C    PHOTON q ----> Q gluon
C
C         P1(q)-----//////--------Q1(q)
C                  ////////
C         P2(PH)-----//////--------Q2(gluon)
C
C
C
C
C for light quarks without masses
c changed for including Q2 of photon
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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEND.
      Double Precision P2E(4),P3E(4),Q1E(4),Q2E(4)
      Double Precision SPHI,PHIT
	Double Precision ALPHAS,ALPHA_S,ALPH_EM,DOT,SH,TH,UH
	Double Precision SUMT,SUML,Y,XP,ZQ,EPSILON,WT1,GF,SIG0t,SIG0s
	Double Precision Sigi,cos2p,cosp,sig1,sig2,sig0,sig
	Integer I,IPHIDEP
      REAL SNGL
      IPHIDEP=1
C.....MARTIX ELEMENT
c      write(6,*) 'here in new eleqcdc '
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
      IF(IPHIDEP.EQ.0) THEN
         IF(INTER.EQ.2) THEN
            write(6,*) ' INTER = 2 not implemented'
            STOP
         ENDIF
         DO 10 I=1,4
            P2E(I) = P(NIA2,I)
            Q1E(I) = P(NF1,I)
            Q2E(I) = P(NF2,I)
   10    CONTINUE
         SH = 2.D0 * DOT(Q1E,Q2E)
         TH = - 2.D0 * DOT(P2E,Q2E)
         UH = - 2.D0 * DOT(P2E,Q1E)

c for transverse photons
         SUMT = 16.D0*8.D0/3.D0*PI*PI*(-TH/SH - SH/TH + 2.D0*UH*Q2/SH/
     +   TH - 2.D0*UH*Q2/(SH+Q2)**2)
c for longitudinal photons
         SUML = 16.D0*32.D0/3.D0*PI*PI*Q2*(Q2+SH+TH)/(SH+Q2)**2
         IF(SUMT.LE.0.0D0.OR.SUML.LT.0.0D0) THEN
            write(6,*) ' eleqcdc ',TH,SH,UH
            write(6,*) ' eleqcdc ',SUMT,SUML
         endif
c         write(6,*) ' eleqcdc SUMT,SUML',SUMT,SUML
c epsilon
         Y = DBLE(YY)
         EPSILON = (1.D0 - Y)/(1.D0 - Y + Y**2 /2)

         WT1= ALPHA_S*ALPH_EM*(SUMT + EPSILON * SUML)
         IF(WT1.LE.0.D0) THEN
            write(6,*) ' eleqcdc ',WT1,EPSILON
            write(6,*) ' eleqcdc ',TH,SH,UH
            write(6,*) ' eleqcdc ',SUMT,SUML
         endif
      ELSEIF(IPHIDEP.EQ.1) THEN
         DO 20 I=1,4
            P2E(I) = p(nia1,I)
            P3E(I) = p(nia2,I)
            Q1E(I) = P(NF1,I)
            Q2E(I) = P(NF2,I)
   20    continue
         y = DBLE(YY)
         SH = 2.D0 * DOT(Q1E,Q2E)
         TH = - 2.D0 * DOT(P3E,Q2E)
         UH = - 2.D0 * DOT(P3E,Q1E)
         xp = Q2/(SH+Q2)
         zq = DOT(P3E,Q1E)/DOT(P3E,P2E)
         zq = 2.D0*DOT(P3E,Q1E)/(SH+Q2)
         WT1=0.d0
c         if(y.ge.1.) write(6,*) ' y = ',y
         SPHI = DLANGL(P(1,1),P(1,2))
         call DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
         phit = DLANGL(P(NF1,1),P(NF1,2))
c      phitt = raphi(1,nia2,nf1)
c     write(6,*) 'eleqqf  phitt,phit ',phittt,phitt,phit
         call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
         cosp = cos(phit)
         cos2p = cos(2.D0*phit)
c fill MEINFO
         ZQGKI = SNGL(zq)
         XPGKI = SNGL(xp)
         if(phit.ge.0.0) then
            PHITGKI = SNGL(PHIT)
         else
            PHITGKI = SNGL(2.D0*PI+PHIT)
         endif
c...........
         sig0t = 2.D0*ALPHA_S/3.D0/pi*((zq**2+xp**2)/
     +    (1.D0-xp)/(1.D0-zq) +
     +   2.D0*(xp*zq+1.D0))
         sig0s = 2.D0*ALPHA_S/3.D0/pi*4.D0*xp*zq
         sigi = 0.D0
         IF(INTER.EQ.2) THEN
            sigi = 2.D0*ALPHA_S/3.D0/pi*((zq**2+xp**2)/ (1.D0-xp)/
     +      (1.D0-zq) + 2.D0*(xp+zq))
         ENDIF

c      write(6,*) ' new sig_t ',sig0t*64*alph_em*pi**3/alph_em/alpha_s
c      write(6,*) ' new sig_l ',sig0s*64*alph_em*pi**3/alph_em/alpha_s

         sig1 = 4.D0*ALPHA_S/3.D0/pi*y*
     +    DSQRT(((1.D0-y)*xp*zq)/(1.D0-xp)/(1.D0-zq))
         sig1 = sig1*((1.D0-2.D0/y)*(1.D0-zq-xp+2.D0*xp*zq))
         sig2 = 4.D0*ALPHA_S/3.D0/pi*((1.D0-y)*xp*zq)
         sig0 = 0.5D0*(1.D0+(1.D0-y)**2)*sig0t + (1.D0-y)*sig0s
     +          +ISIGN(1,K(1,2))*y*(1.d0-0.5d0*y)*sigi
         sig = sig0+ cosp * sig1 + cos2p * sig2
c include 64/pi*alph*pi**3 to compare with s_hat expression and
c integration over angles
         IF(INTER.LT.2) THEN
            wt1 = alph_em**2 /pi/Q2/y * sig *64.D0*pi**3
         ELSEIF(INTER.EQ.2) THEN
            GF = PI*ALPH_EM/(SIN2W*XMW2*DSQRT(2.D0))
ccc            wt1 = GF**2/pi/Q2/y * sig *64.D0*pi**3
            wt1 = GF**2 *Q2/y * 8.D0 * sig
c is that really correct ????? to be checked
            wt1 = wt1 /(1.D0 + Q2/XMW2)**2
         ENDIF
c      write(6,*) 'eleqcdc:',zq,xp,yy,sig,alph_em,pi
      ENDIF
      RETURN
      END
