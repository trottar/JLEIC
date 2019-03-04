*CMZ :  2.08/00 06/06/99  15.53.20  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  18.47.29  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE ELEQQF(WT1)
C
C    PHOTON GLUON ----> Q Q_BAR
C
C         P1(G)-----//////--------Q1(q)
C                  ////////
C         P2(PH)-----//////--------Q2(q_bar)
C
C    full martix element G.Schuler
C    e g --> e' q q_bar
C
C    Notation: p = P3E
C              q = P2E
C              l = P1E
C              pf = Q1E
C              pf' = Q2E
C              sh = s_hat = gamma glu CM energy
C              sg = e glu CM energy
C              x = Q2/(2P.q)
C              xg = Q2/(2p.q)
C              y = p.q/p.l
C              z = p.pf/p.q
C              |M|**2 = g**2 * eq**2/Q**4 *sum gf L_munu *H^munu
      IMPLICIT NONE
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


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
      Double Precision P0E(4),P2E(4),P3E(4),Q1E(4)
      Double Precision G1(2),G2(2),GA(2),GB(2)
      Double Precision gff(4)
	Double Precision Vff,QF2
      Double Precision SPHI,phit,mf,mfp
      Double Precision ksiw,raphi,wt1,y,sg,sum,xg,m2plus,m2minus,m0
	Double Precision z,cosp,cos2p,f,g
	Double Precision H11,H12,H21,H22,H41,H14,H24,H34,H61
	Double Precision G3,GC,deltal,rho
	Double Precision ALPHA_S,ALPH_EM,DOT,ALPHAS
	Integer Luchge,I
	Double Precision Wmax
	Integer IMIX,IGENFL
      COMMON /OALPHAS/ WMAX,IMIX
      COMMON/GENWEI/IGENFL
      REAL SNGL
      external raphi
      IF(INTER.LT.2) THEN
c charge for charm
         IF(IPRO.EQ.14) THEN
            IF(AM(1).LT.2.D0) QF2 = 4.D0/9.D0
            IF(AM(1).GE.2.D0) QF2 = 1.D0/9.D0
c         write(6,*) ' ELEQQF IPRO=14',QF2
         ELSEIF(IGENFL.EQ.1) THEN
            QF2 = DFLOAT(LUCHGE(KPA))**2/9.D0
c         write(6,*) ' ELEQQF IGENFL=1 KPA ',KPA,QF2

         ELSE
c sum of charges for u d s = (2/3)**2 + (1/3)**2 + (1/3)**2 = 6/9 = 2/3
            QF2=2.D0/3.D0
c         write(6,*) ' ELEQQF SONST KPA ',KPA,QF2
         ENDIF
      ELSEIF(INTER.EQ.2) THEN
c charge for charm
         IF(IPRO.EQ.14) THEN
            Vff = 1.D0
         ELSEIF(IGENFL.EQ.1) THEN
            Vff = 1.D0
         ELSE
c sum of charges for u d s = 3
            Vff=3.D0
c         write(6,*) ' ELEQQF SONST KPA ',KPA,QF2
         ENDIF
      ENDIF
      WT1 = 0.0D0
   10 CONTINUE
      if(IPRO.EQ.13) THEN
         AM(1) = 0.D0
         AM(2) = 0.D0
      ENDIF
      mf = AM(1)
      mfp = AM(2)
c         write(6,*) 'masses ',mf,mfp,ipro
c
c      IF(IPRO.EQ.14) THEN
c         write(6,*) 'masses ',mf,mfp
c         write(6,*) ' charge ',DSQRT(QF2)
c         ENDIF
      DO 20 I=1,4
         P0E(I) = P(1,I)
         P2E(I) = p(nia1,I)
         P3E(I) = p(nia2,i)
         Q1E(I) = P(NF1,I)
   20 continue
      y = DBLE(YY)
c      sgp = DBLE(XPR)*SSS
c check for pomeron_gluon
      sg= 2.D0*DOT(P0E,P3E)
c      call DULIST(1)
c      xgp = Q2/y/sgp
      XG = Q2/2.D0/DOT(P2E,P3E)
c      write(6,*) ' eleqqf: sg",sg,xg",xg ',sgp,sg,xgp,xg
      m2plus = (mf**2 + mfp**2)/y/sg
      m2minus = (mf**2 - mfp**2)/y/sg
      m0      = mf*mfp/y/sg
      z = DOT(P3E,Q1E)/DOT(P3E,P2E)
ctest
c      ZQMIN = 0.5D0 - DSQRT(0.25D0 - PT2CUT(IPRO)*XG/Q2/(1.D0-XG))
c      write(6,*) ' eleqqf: zqmin,z ',zqmin,z,pt2cut(ipro)
c      SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
c      CHIMIN=(SMIN+Q2)/(Y*SSS)
c      CHIMIN = DMAX1(CHIMIN,Q2/Y/SSS)
c      x = q2/y/sss
c      XPMIN = X
c      XPMAX = X/CHIMIN
c      write(6,*) y,sss,q2,chimin,x
c      write(6,*) 'eleqqf xpmax,xpmin,xp ',xpmax,xpmin,xg
ctest
c      write(6,*) ' eleqqf z',z,nia1,nia2,nf1,nf2
ctest      if(y.ge.1.) write(6,*) ' y = ',y
      SPHI = DLANGL(P(1,1),P(1,2))
c     if(igenfl.eq.1) then
c      call dulist(1)
c      endif
      call DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      phit =  DLANGL(P(NF1,1),P(NF1,2))
c     if(igenfl.eq.1) then
c      call dulist(1)
c      write(6,*) ' phi ',phit
c      endif
c      phitt = raphi(1,nia2,nf1)
c     write(6,*) 'eleqqf  phit ',phit
      call DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
c we are in a system with g(p) in - z and gamma in + z
c therefore -cos(phit)
      cosp = DBLE(-cos(phit))
      cos2p = DBLE(cos(2.D0*phit))
c fill MEINFO
      XPGKI = sngl(xg)
      ZQGKI = sngl(z)
      if(phit.ge.0.0) then
         PHITGKI = SNGL(PHIT)
      else
         PHITGKI = SNGL(2.D0*PI+PHIT)
      endif

c............
      f = (1.D0-xg)*z*(1.D0-z) + (z*(mf**2 - mfp**2) - mf**2)/y/sg
      g = (1.D0-xg)*(1.D0-z) + xg*z + m2minus
      if(f*xg.le.0.) then
c         write(6,*) 'eleqqf f,xg,f*xg,z',f,xg,f*xg,z,mf**2/y/sg
c         write(6,*) 'eleqqf (1.D0-xg)*z*(1.D0-z)',(1.D0-xg)*z*(1.D0-z)
c         write(6,*) 'eleqqf mf**2/y/sg',mf**2/y/sg
         return
      endif
      H11 = 4.D0/z/(1.D0-z)*(2.D0*z*(1.D0-z) + 2.D0*xg*(1.D0-xg) -1.D0)
     +    + 4.D0/z**2/(1.D0-z)**2 *(-m2plus**2 + m2plus*m2minus*
     +      (2.D0*z-1.D0) + m2plus*(2.D0*z*(z-1.D0)*(xg-1.D0) -xg) +
     +        m2minus*xg*(2.D0*z-1.D0))
      H21 = 16.D0*(xg/z/(1.D0-z) + (m2plus +
     + (1.D0-z)*m2minus)/z**2/(1D0-z))
      H41 = 32D0*xg/z/(1.D0-z) + 16.D0/z**2/(1.D0-z)**2 *
     +      (m2plus + (1.D0-2.D0*z)*m2minus)
      H61 = - 32.D0*xg/z/(1.D0-z) - 32.D0/z**2/(1.D0-z) *
     +      (m2plus + (1.D0-z)*m2minus)
      H12 = - 16D0* m0 *(1.D0-xg)/z/(1.D0-z) +
     +       8.D0*m0/z**2/(1.D0-z)**2 *
     +      (m2plus + (1.D0-2.D0*z)*m2minus)
      H22 = - 32.D0*m0/z/(1.D0-z)
      H14 = - 8.D0/z**2/(1.D0-z)*(-z**2+z*xg+m2plus+(1.D0-z)*m2minus)
      H24 = 0.D0
      H34 = 8.D0*(1.D0-2.D0*xg)/z/(1.D0-z) -
     +      8.D0*(m2plus+(1.D0-2.D0*z)*m2minus)/z**2/(1.D0-z)**2
      G1(1) = -2.D0*H11 + f*H41
      G2(1) = 2.D0*xg*f*H41 + H21 + g*(H61 + g * H41)
      GA(1) = DSQRT(xg*f) * (H61 + 2.D0*g*H41)
      GB(1) = f * H41
      G1(2) = - 2.D0* H12
      G2(2) = H22
      GA(2) = 0.d0
      GB(2) = 0.d0
      G3 = H14 + z*H24 -g*H34
      GC = DSQRT(xg*f)*(H24 - 2.D0*xg*H34)
      IF(INTER.EQ.0) THEN
         gff(1) = QF2
         gff(2) = QF2
         gff(4) = 0.D0
         deltal = 0.5d0
      ELSEIF(INTER.EQ.2) THEN
         ksiw = 1.D0/(SIN2W*8.D0)*Q2**2/XMW2**2/(1.D0 + Q2/XMW2)**2
         rho = 0.d0
         gff(1) = 4.D0*Vff**2*ksiw**2*(1.D0-ISIGN(1,K(1,2))*rho)
         gff(2) = 0.d0
         gff(4) = 4.D0*Vff**2*ksiw**2*(ISIGN(1,K(1,2)) - rho)
         deltal = 1.d0
      ENDIF
      SUM = gff(1)*(xg*y**2*G1(1) + (1.D0 - y)*G2(1) +(2.D0-y)*
     +DSQRT(1.D0-y)* cosp*GA(1) + 2.D0*xg*(1.D0-y)*cos2p*GB(1)) +gff(2)
     +*(xg*y**2*G1(2) + (1.D0 - y)*G2(2)) +gff(4)*(xg*y*(2.D0-y)*G3+
     +2.D0*y*DSQRT(1.d0-y)*cosp*GC)

C.....MARTIX ELEMENT
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
C this is for z,phi integrated and full differential
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
c       write(6,*) 'elleqqf ',alpha_s,SUM,qf2,y,q2,alph_em,irunaem
c      write(6,*) ' sum,alpha_s,q2q ',sum,alpha_s,q2q
      SUM = SUM *ALPHA_S* deltal*2.D0*ALPH_EM**2
      SUM = SUM / 16.D0 / y /Q2**2
c this is correction for 2 body phase space
      SUM = SUM*2.D0/pi
      WT1=SUM
c      write(6,*) WT1,ksiw
      RETURN
      END
