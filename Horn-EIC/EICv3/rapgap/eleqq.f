*CMZ :  2.08/00 06/06/99  17.46.35  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  17.27.49  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE ELEQQ(WT1)
C
C    e p  ----> Q Q_BAR  p
C
C
      implicit none
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
      DOUBLE PRECISION PK1(4),PK2(4)
      DOUBLE PRECISION DBW2(4),DBMX(4)
      DOUBLE PRECISION STHETA,SPHI,sphi1
      DOUBLE PRECISION MX2,M2
      DOUBLE PRECISION WMAX
      DOUBLE PRECISION I_L,I_T
      INTEGER IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      INTEGER IGENFL
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION XPOM
      COMMON/BARTELS/XPOM
      REAL SNGL
      DOUBLE PRECISION DOT
      INTEGER I
      DOUBLE PRECISION   sst,betat,xppt,t,qqt
      common   /parameter/ sst,qqt,betat,xppt,t
      DOUBLE PRECISION WT1,sigdt,sigdl,sigdi,sigda
      DOUBLE PRECISION ALPH_EM,Y,W2,XB,PT2
      DOUBLE PRECISION PTFACT,mf,ksi,ALPHA_S,ALPHAS
      DOUBLE PRECISION pt2max,charge,const
      DOUBLE PRECISION phi
      DOUBLE PRECISION cosp,cos2p,raphi

      common   /phi/       phi
      INTEGER IPHI
      COMMON /SEL/IPHI
      DOUBLE PRECISION D_XGX,glu,dglu,xglu
      EXTERNAL ALPHAS,D_XGX,XGLU,raphi
      t = dble(t2gki)
      sst = SSS
      qqt = q2
      xppt = xpom
      iphi=0
      WT1 = 0.D0
      y = DBLE(YY)
      DBMX(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBMX(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBMX(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBMX(4)=  P(NIA1,4) + P(NIA1+1,4)
      MX2=DOT(DBMX,DBMX)
c  goto gamma p system
      DBW2(1)=  P(NIA1,1) + P(2,1)
      DBW2(2)=  P(NIA1,2) + P(2,2)
      DBW2(3)=  P(NIA1,3) + P(2,3)
      DBW2(4)=  P(NIA1,4) + P(2,4)
      W2=DOT(DBW2,DBW2)
      xb = Q2/(W2 + Q2)

      betat = xb/xpom

c  goto gamma pom system
      DBW2(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBW2(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBW2(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBW2(4)=  P(NIA1,4) + P(NIA1+1,4)
c     write(6,*) ' before boost '
c      call dulist(1)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  -DBW2(1)/DBW2(4),-DBW2(2)/DBW2(4),
     +  -DBW2(3)/DBw2(4))
c     write(6,*) ' before phi '
c      call dulist(1)
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
c     write(6,*) ' before theta '
c      call dulist(1)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c     write(6,*) ' after boost '
c      call dulist(1)

      DO 10 I=1,4
         PK2(I) = P(NF1,I)
         PK1(I) = P(NF2,I)
   10 continue
      SPHI1 = DLANGL(P(1,1),P(1,2))
      call DUDBRB(0,0,0.D0,-sphi1,0.d0,0.d0,0.d0)
      phi =  DLANGL(P(NF1,1),P(NF1,2))
c      phitt = raphi(1,nia2,nf1)
c     write(6,*) 'eleqqf  phitt,phit ',phittt,phitt,phit
      call DUDBRB(0,0,0.D0,sphi1,0.d0,0.d0,0.d0)
      if(phi.ge.0.0) then
         PHITGKI = SNGL(PHI)
      else
         PHITGKI = SNGL(2.D0*PI+PHI)
      endif

      cosp = dcos(phi)
      cos2p = dcos(2.d0*phi)
      m2 = DOT(PK2,PK2) + 2.D0*DOT(PK2,PK1) + DOT(PK1,PK1)
      m2 = mx2
      PT2 = DBLE(PT2H)
      y = DBLE(YY)
      IF(p(nf1,5).lt.1.d0) then
         mf = 0.0d0
         charge = 2.d0/3.d0
      elseif(p(nf1,5).gt.1.d0.and.p(nf1,5).lt.4.d0) then
         mf = p(nf1,5)
         charge = 4.d0/9.d0
      endif
      ksi = (Q2 + mf**2/pt2*(Q2+m2))/m2
      glu = xglu(xpom,pt2*(1.d0+ksi))
      dglu = 0.d0
      dglu = d_xgx(xpom,pt2*(1.d0+ksi))
calculate I_T
      I_T = 2.d0/(1.d0 + ksi)**3 *(
     +      2.d0*ksi/pt2  * glu +
     +      (1.d0 -ksi -2.d0*ksi*dlog(ksi/(1.d0+ksi)))*
     +       dglu)
calculate I_L
      I_L = (pt2 + mf**2)*Q2/pt2/m2/(1.d0 + ksi)**3
     +      *((ksi - 1.d0)/pt2 * glu +
     +        (2.d0 +(1.d0-ksi)*dlog(ksi/(1.d0+ksi)))*
     +         dglu)
      if((PT2+mf**2)/m2.GT.0.25D0) GOTO 20
c jacobian dM**2 --> d x_pom
c now in partdh included
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
      Q2Q = pt2*(1.d0+ksi)
      ALPHA_S=ALPHAS(DSQRT(Q2Q))
      pt2max=m2/4.d0 - 4.d0*DBLE(ULMASS(3))**2
      if(pt2.ge.pt2max) goto 20
      PTFACT = (PT2 + mf**2)/DSQRT(1.D0 - 4.d0*(PT2+mf**2)/m2)
      const = alph_em*pi**2*alpha_s**2*charge

      sigdt = 1.d0/m2**2/pt2/12.d0*const*
     +  PTFACT*( (1.d0 - 2.d0*(pt2+mf**2)/m2)*I_T**2
     + + mf**2*pt2*m2**2/(pt2+mf**2)/Q2**2 * I_L**2
     + )
      sigdl = 1.d0/m2**2/q2*4.d0/3.d0*ptfact*const*I_L**2
      sigdi = 1.d0/m2**2/dsqrt(pt2*q2)/3.d0*const*(PT2 + mf**2)*I_T*I_L
      sigda = 1.d0/m2**3/pt2/12.d0*ptfact*const*(PT2 + mf**2)*I_T**2
c     write(6,*) ' sigs ',sigdt,sigdl,sigdi,sigda
c     write(6,*) ' I_T,I_L ',I_T,I_L,ksi
      wt1 = alph_em/2.d0/y/Q2/pi**2*(
     +      (1.d0 + (1.d0 - y)**2)/2.d0 * sigdt
     +       + (1.d0 - y) * sigdl - 2.d0*(1.d0-y)*cos2p*sigda
     +       + (2.d0 -y)*dsqrt(1.d0 - y)*cosp*sigdi)
c include t dependence a la DL
c now in partdh included
      IF(WT1.LT.0) THEN
         write(6,*) ' wt1  ',wt1
      ENDIF
      IF(WT1.NE.WT1) THEN
         write(6,*) ' wt1  ',wt1
         write(6,*) ' sigblw:sigdt,sigdl',sigdt,sigdl
         write(6,*) sigdt,sigdl,sigdi,sigda,ptfact,const,m2,I_T,I_L,
     +   PT2
         call dulist(1)
      ENDIF
      CALL DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      RETURN
   20 WT1 = 0.D0
c      write(6,*) ' eleqq: set wt1=0'
      CALL DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      RETURN
      END
