*CMZ :  2.08/04 17/01/2000  13.14.18  by  Hannes Jung
*CMZ :  2.08/02 08/09/99  17.06.07  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.57  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  17.26.32  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE ELEQQG(WT1)
C
C    e p  ----> Q Q_BAR GLUON p
C
C         P1(G)-----//////--------Q1(q)
C                  ////////
C         P2(PH)-----//////--------Q2(q_bar)
C
C    full martix element
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

*KEEP,RGPQCDPOM.
      Integer Iqqg
	COMMON/pqcdpom/ Iqqg
C     SAVE


*KEND.
      DOUBLE PRECISION P0E(4),P2E(4),P3E(4)
      DOUBLE PRECISION PK1,PK2,PK3,PL
      COMMON /QQG/ PK1(4),PK2(4),PK3(4),PL(4)
      DOUBLE PRECISION QPR(4),PP(4)
      DOUBLE PRECISION a1ma1,Q2_c
      COMMON /QQG_CON/a1ma1,Q2_c
      DOUBLE PRECISION DBW2(4),DBMX(4)
      DOUBLE PRECISION STHETA,SPHI
      DOUBLE PRECISION MX2,M2,K12,K1K2,K22
      DOUBLE PRECISION WMAX
      INTEGER IMIX
      COMMON /OALPHAS/ WMAX,IMIX
      INTEGER IGENFL
      COMMON/GENWEI/IGENFL
      DOUBLE PRECISION XPOM
      COMMON/BARTELS/XPOM
      DOUBLE PRECISION WT1,sigdt_p,sigdt_m,sigdl,sigdi,WT2
      DOUBLE PRECISION alpha1,beta2,PREF,ALPH_EM,ALP,Y,W2,XB
      DOUBLE PRECISION m2bs,alpha2,beta1,a21ma2
      REAL SNGL
      REAL SLO,SHI,XEPS,EPS
      REAL SUMM11,XM11,SUMM21,XM21
      REAL SUMM12,XM12,SUMM22,XM22
      REAL SUMM1,XM1
      REAL SUMM2,XM2
      REAL FIL,FIU
      REAL L2MIN,L2MAX
      DOUBLE PRECISION DOT
      INTEGER I,iphase
      COMMON /QQG_C/L2MIN,L2MAX
      Double Precision   sst,betat,xppt,t,qqt
      common   /parameter/ sst,qqt,betat,xppt,t
      Double Precision DK,th,uh,qf2,k22fr
      EXTERNAL SUMM11,SUMM21,SUMM12,SUMM22,SUMM1,SUMM2
      EXTERNAL FIL,FIU,DOT,DK
      Double Precision D_XGX,glu,xglu,ksi,xcorr,alphas
      Double Precision FACQ2,beta,z,arg,Coeff,m2cc,mq2
      EXTERNAL ALPHAS,D_XGX,XGLU
cc iphase = 1 only phase space
c      iphase = 1
c iphase = 0 x section
      iphase = 0
      t = dble(t2gki)
      WT1 = 0.D0
      WT2 = 0.d0
      XEPS = 0.05
      SLO = 0.
      SHI = 1.
      L2MIN = 1.
      L2MIN = 0.3
      L2MAX = 1000.
      Q2_c = Q2
      WT1 = 0.0D0
      y = DBLE(YY)
      Xb = Q2/Y/SSS
c         xb = Q2/(W2 + Q2)
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
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
c  goto gamma pom system
      DBW2(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBW2(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBW2(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBW2(4)=  P(NIA1,4) + P(NIA1+1,4)

      CALL DUDBRB(0,0,0.D0,0.D0,
     +  -DBW2(1)/DBW2(4),-DBW2(2)/DBW2(4),
     +  -DBW2(3)/DBw2(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c      call dulist(1)
      DO 10 I=1,4
         P0E(I) = P(1,I)
         P2E(I) = p(nia1,I)
         P3E(I) = p(nia2,I)
         PK3(I) = P(NF1,I)
         PK2(I) = P(NF1+1,I)
         PK1(I) = P(NF2,I)
         QPR(I) = P(NIA1,I) + XB*P(2,I)
         PP(I) = P(2,I)
   10 continue
      th = -Q2 - 2.d0*DOT(P2E,PK3) + P(NF1,5)**2
      uh = -Q2 - 2.D0*DOT(P2E,PK1) + P(NF2,5)**2
      beta2 = 2.D0*DOT(PK2,QPR)/(W2+Q2)
      beta1 = 2.D0*DOT(PK1,QPR)/(W2+Q2)
c     write(6,*) 'q"2 ',dot(qpr,qpr)
      alpha1 = 2.D0*DOT(PK1,PP)/(W2+Q2)
      alpha2 = 2.D0*DOT(PK2,PP)/(W2+Q2)
      m2 = DOT(PK3,PK3) + 2.D0*DOT(PK3,PK1) + DOT(PK1,PK1)
      m2 = 2.D0*DOT(PK3,PK1)
c      call dulist(1)
      if(m2.lt.1.d0) goto 50
      K12 = PK1(1)**2 + PK1(2)**2
      K22 = PK2(1)**2 + PK2(2)**2
      IF(IHFLA.GE.4) THEN
         K12 = K12 + P(NF2,5)**2
      ENDIF

      K1K2 = (PK1(1)+PK2(1))**2 + (PK1(2) + PK2(2))**2


      m2bs = MX2 - beta2 * W2
c alp = alphas fixed
      ALP = 0.25D0

      IF(alpha1.gt.1.D0.or.alpha1.lt.0.D0) THEN
         write(6,*) ' alpha1 = ',alpha1
c      write(6,*) ' K1**2 = ',K12
c      write(6,*) ' (K1+K2)**2 = ',K1K2
c      write(6,*) ' M_X**2, m**2,W**2 ',MX2,m2,W2
c      call dulist(1)
         GOTO 40
      ENDIF
c a1ma1 = alpha(1-alpha)
      a1ma1 = alpha1*(1.d0 -alpha1)
c a21ma2 = alpha^2 + (1-alpha)^2
      a21ma2=alpha1**2 + (1.d0-alpha1)**2
c now here the rest of cuts...
      If(Iqqg.eq.0) Then
         if(dabs(th).lt.pt2cut(ipro)) goto 40
         if(dabs(uh).lt.pt2cut(ipro)) goto 40
      Endif
c this cut is essential for the approximation used
      if(pk2(3).GT.0.D0) goto 40
c
      if (iphase.eq.1) goto 30
c select here approximate formula or full integration
      If(Iqqg.ge.1) then
         XM11 = 0.
         Xm22 = 0.
         Xm12 = 0.
         Xm21 = 0.
         XM1 = 0.
         XM2 = 0.
c cuts on pt
         if(K12.lt.pt2cut(ipro).AND.IHFLA.LE.3) goto 20
         if(k12.GT.m2/4.d0) goto 20
         if(k22.lt.0.001) goto 20
c alp = alphas running
         Q2Q = K12
         ALP = ALPHAS(DSQRT(Q2Q))
         k22fr = k22+0.04D0
         m2 = m2bs
cccc         ksi = dmax1(k22fr,1.d0)
         ksi = dmax1(k22fr*(MX2+Q2)/(MX2-m2),1.d0)
cccc         ALP = ALPHAS(DSQRT(ksi))
         ALP = ALPHAS(DSQRT(k12))
         glu = xglu(xpom,ksi)
         Facq2 = (Q2**2+m2**2)/(Q2+m2)**4
c        write(6,*) ' stand. Facq2 ',Facq2
         mq2 = 0.d0
         IF(IHFLA.ge.4) then
            z= (Q2 + m2)/(Q2 + MX2)
            beta = Q2/(Q2 + MX2)
            arg = beta/z
            Facq2= 2.d0*arg**2/Q2**2
c test with standard AP splitting fct
            Coeff = (arg**2 + (1.d0 -arg)**2)/2.d0
c now with heavy quark splitting fct
c acc. eq.(16) in Ellis,Sexton Nucl.Phys.B282(1987) 642
            m2cc = DOT(PK3,PK3) + 2.D0*DOT(PK3,PK1) + DOT(PK1,PK1)
            Coeff = (arg**2 + (1.d0 -arg)**2 + 2.d0*AM(1)**2/m2cc)/
     +      2.d0
c            rr=AM(1)**2/Q2
c        sqrtv=dsqrt(1.d0-4.d0*am(1)**2*arg/Q2/(1.d0-arg))
c        Coeff = (arg**2 + (1.d0 -arg)**2 +
c     &  4.d0*arg*(1.d0-3.d0*arg)*rr )*
c     &  dlog((1.d0+sqrtv)/(1.d0-sqrtv)) +
c     &      sqrtv*(-1.d0 + 4.d0*arg*(1.d0 - arg) -
c     &      4.d0*arg*(1.d0-arg)*rr)
            Facq2 = Facq2 * Coeff
            ALP = ALPHAS(DSQRT(k12)+2.d0*AM(1))
            mq2 = AM(1)**2
c        write(6,*) mq2,am(1)
c           write(6,*) ' new   Facq2 ',Facq2
         ENDIF
         XM11= SNGL(8.D0* Facq2 * (m2/k12+mq2)**2 /k22fr**2)
c add additional term
         IF(Iqqg.eq.2) Then
            XCORR = ((MX2 - m2)/(MX2+Q2))**4 * ((MX2+2.d0*m2+3.d0*Q2)/
     +      (MX2+Q2))**2
         ELSE
            XCORR = 1.D0
         ENDIF
         XM11 = XM11*SNGL(glu**2 * XCORR)
c apply a factor of 2 because we have two solutions to eq.(2.18)
         XM11 = 2.*XM11
         XM11 = SQRT(XM11)
         XM22 = XM11
c the following lines are essential for the approximation used here
c a1ma1 = alpha(1-alpha)
         a1ma1 = (k12+mq2)/m2
c a21ma2 = alpha^2 + (1-alpha)^2
         a21ma2 = 1.D0 - 2.D0*(k12+mq2)/m2
   20    Continue
      elseif(iqqg.eq.0) then
c calculate integral over l
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM11,EPS,XM11)
c      write(6,*) ' integral XM11 ',XM11
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM21,EPS,XM21)
c      write(6,*) ' integral XM21 ',XM21
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM12,EPS,XM12)
c      write(6,*) ' integral XM12 ',XM12
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM22,EPS,XM22)
c      write(6,*) ' integral XM22 ',XM22
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM1,EPS,XM1)
c      write(6,*) ' integral XM1 ',XM1
         EPS = XEPS
         CALL GADAP2(SLO,SHI,FIL,FIU,SUMM2,EPS,XM2)
c      write(6,*) ' integral XM2 ',XM2
      endif
   30 continue
c new prefactor for the use with lorentz invariant phase space
c according to mark
      Pref = 9.d0/32.D0/PI
c sum over light charges = 2/3
      IF(IHFLA.LE.3) THEN
         QF2 = 2.D0/3.D0
      ELSEIF(IHFLA.EQ.4) THEN
         QF2 = 4.D0/9.D0
      ELSEIF(IHFLA.EQ.5) THEN
         QF2 = 1.D0/9.D0
      ENDIF
      Pref = Pref*QF2*ALPH_EM*ALP**3
      sigdt_p = Pref*a21ma2*a1ma1* dble(XM11*
     +XM11+XM12*XM12+XM21*XM21+XM22*XM22)

      sigdt_m = Pref*a1ma1**2*
     +dble(XM11*XM11+XM12*XM12-XM21*XM21-XM22*XM22)
      sigdl = Pref*4.D0*a1ma1**3*Q2*dble(XM1*XM1+XM2*XM2)
      sigdi = Pref*a1ma1**2*(1.D0-2.D0*alpha1)*0.5D0*DSQRT(Q2)*
     +DBLE(XM11*XM1+XM12*XM2+XM1*XM11+XM2*XM22)
      if(iphase.eq.1) then
         sigdt_p =1.D0
         sigdt_m =0
         sigdl = 0
         sigdi = 0
      endif
c      write(6,*) ' pi,m2,beta2 ',pi,m2,beta2
c      write(6,*) ' Pref ',Pref,a1ma1
c      write(6,*) 'sigdt_p ',sigdt_p
c      write(6,*) 'sigdt_m ',sigdt_m
c      write(6,*) 'sigdl ',sigdl
c      write(6,*) 'sigdi ',sigdi
      wt1 = alph_em/2.d0/y/Q2/pi**2*( (1.d0 + (1.d0 - y)**2)/2.d0 *
     +sigdt_p - 2.d0*(1.d0 - y) * sigdt_m + (1.d0 - y) * sigdl + (2.d0
     +-y)*dsqrt(1.d0 - y)*sigdi)
c         write(6,*) ' wt1  ',wt1
      IF(WT1.LT.0) THEN
	   wt1  = 0.
c         write(6,*) ' wt1  ',wt1
c         write(6,*) ' a21ma2 ',a21ma2,a1ma1
c         write(6,*) 'sigdt_p ',sigdt_p
c         write(6,*) 'sigdt_m ',sigdt_m
c         write(6,*) 'sigdl ',sigdl
c         write(6,*) 'sigdi ',sigdi
      ENDIF
      IF(WT1.ne.wt1) THEN
         write(6,*) ' error wt1  ',wt1
      ENDIF
      WT2 = WT2 + WT1
   40 CONTINUE
   50 CONTINUE
c      write(6,*) ' wt1  ',wt1
c      write(6,*) ' end eleqqg '
c      call dulist(1)
      CALL DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBW2(1)/DBW2(4),DBW2(2)/DBW2(4),
     +  DBW2(3)/DBW2(4))
      WT1 = WT2
      RETURN
      END
