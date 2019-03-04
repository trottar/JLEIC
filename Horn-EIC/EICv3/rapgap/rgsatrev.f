*CMZ :  2.08/05 27/03/2000  15.59.09  by  Hannes Jung
*CMZ :  2.08/04 22/12/99  15.39.29  by  Hannes Jung
*CMZ :  2.08/02 02/11/99  12.47.44  by  Hannes Jung
*CMZ :  2.08/00 14/06/99  19.13.58  by  Hannes Jung
*-- Author :    Hannes Jung   30/05/99
      Subroutine RGSATREV
      Implicit None

      Double Precision  LIMKIND
      Integer           IKIND,      Nevent, Ievent
      Common /CKIND/    LIMKIND(20),IKIND , Nevent, Ievent
      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC

      Double Precision    DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT
      Double Precision    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT
      COMMON /CCOMPCROS/  DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT,
     &                    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT

      Double Precision  EPS
      Common  /mysatrap/EPS
      Double Precision BOCHCK,SPHI,STHETA,PEP,PEZ,PEG,PEGZ
      Double Precision POMDGA,EN,PZC,PT,t2,CPHP,SPHP,spom,PHIP
      Integer NPFIN,NDFF,KI
      Double Precision RN,DOT
      Double Precision  Pi
      Double Precision draprn
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

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


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

      PI = 3.1415927
      YY = Ybj
      XEL = Ybj
      XPR = Xbj
      T2GKI = Tdf
      XFGKI = xpom

      DLTOT = DLQQT+DLQQGT + EPS*(DLQQL+DLQQGL)
      DLCTOT= DCQQT+DCQQGT + EPS*(DCQQL+DCQQGL) + DLTOT

      LIMKIND(1) = DLQQT/DLCTOT
      LIMKIND(2) = LIMKIND(1) + EPS*DLQQL/DLCTOT
      LIMKIND(3) = LIMKIND(2) +     DLQQGT/DLCTOT
      LIMKIND(4) = LIMKIND(3) + EPS*DLQQGL/DLCTOT

      LIMKIND(5) = LIMKIND(4) +      DCQQT/DLCTOT
      LIMKIND(6) = LIMKIND(5) +  EPS*DCQQL/DLCTOT
      LIMKIND(7) = LIMKIND(6) +      DCQQGT/DLCTOT
      LIMKIND(8) = LIMKIND(7) +  EPS*DCQQGL/DLCTOT


      IKIND = 20
      RN = draprn()
      If(RN.le.LIMKIND(1))                      IKIND = 1
      If(RN.gt.LIMKIND(1).and.RN.le.LIMKIND(2)) IKIND = 2
      If(RN.gt.LIMKIND(2).and.RN.le.LIMKIND(3)) IKIND = 3
      If(RN.gt.LIMKIND(3).and.RN.le.LIMKIND(4)) IKIND = 4

      If(RN.gt.LIMKIND(4).and.RN.le.LIMKIND(5)) IKIND = 11
      If(RN.gt.LIMKIND(5).and.RN.le.LIMKIND(6)) IKIND = 12
      If(RN.gt.LIMKIND(6).and.RN.le.LIMKIND(7)) IKIND = 13
      If(RN.gt.LIMKIND(7).and.RN.le.LIMKIND(8)) IKIND = 14

      MQUARK = 0.0001
      If(ikind.gt.10) Then
         MQUARK = MCHARM
      Endif
      If(IKIND.eq.20) Then
         Write(*,*) 'Error IKIND '
      Endif
c     write(6,*) ' dulist '
c     write(6,*) ' in gamma p '
c     call dulist(1)
      NDFF = 4
      NF1 = NIA1+2
      NF2 = NIA1+3
      if(ikind.eq.3.or.ikind.eq.4.or.ikind.eq.13.or.ikind.eq.14) Then
         NF2 = NIA1+4
         Ndff = 5
      Endif
c     write(6,*) ' ikind ',ikind,' ndff ',ndff
      NPFIN=NIA1+NDFF
      N=NPFIN
      K(NPFIN,1)=1
      K(NPFIN,2)=KP
      IF(IABS(KINT(2,2)).EQ.211) K(NPFIN,2) = 2112
      K(NPFIN,3)=2
      P(NPFIN,5)=DBLE(ULMASS(K(NPFIN,2)))
c for construction of  pomeron goto gamma p system
      DBCMS(1)=  P(NIA1,1) + P(2,1)
      DBCMS(2)=  P(NIA1,2) + P(2,2)
      DBCMS(3)=  P(NIA1,3) + P(2,3)
      DBCMS(4)=  P(NIA1,4) + P(2,4)
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +              + (DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 20
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  -DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +  -DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(0,0,-STHETA,0.D0,0.d0,0.d0,0.d0)
c new
      PEP = DSQRT(P(2,3)**2 + DBLE(ULMASS(KP)**2))
      P(2,4) = PEP
      PEZ = ABS(P(2,3))
      PEG = P(NIA1,4)
      PEGZ = DABS(P(NIA1,3))
c     write(6,*) ' y = ',ybj,' xpom ',xpom
cres      POMDGA = XR * YX * SSS/2.D0
      POMDGA = xpom * Ybj * SSS/2.D0
      T2 = -Tdf
      EN= PEZ*(PEGZ*PEZ+PEG*PEP)/(PEZ*PEG+PEP*PEGZ) - ((T2-P(2,5)**2-
     +P(NPFIN,5)**2)*PEGZ/2.D0 + POMDGA*PEZ)/(PEZ*PEG+PEP*PEGZ)
      P(NPFIN,4) = EN
      PZC =       - (PEP*(PEGZ*PEZ+PEG*PEP)/(PEZ*PEG+PEP*PEGZ) +
     +             ((T2-P(2,5)**2-P(NPFIN,5)**2)*PEG/2.D0
     +             - POMDGA*PEP)/(PEZ*PEG+PEP*PEGZ))
      P(NPFIN,3) = PZC
      PT = (-(T2-P(NPFIN,5)**2)*(PEGZ*PEZ+PEG*PEP)**2 + (T2-P(2,5)**2-
     +P(NPFIN,5)**2)**2*Q2/4.D0 -P(NPFIN,5)**2*POMDGA**2 +(T2+P(2,5)**
     +2-P(NPFIN,5)**2)*POMDGA *(PEGZ*PEZ+PEG*PEP) )/(PEZ*PEG+PEGZ*PEP)*
     +*2 - P(2,5)**2
      PHIP=2.D0*PI*draprn()
      SPHP=DSIN(PHIP)
      CPHP=DCOS(PHIP)
      P(NPFIN,1) = DSQRT(DMAX1(0.D0,PT))*CPHP
      P(NPFIN,2) = DSQRT(DMAX1(0.D0,PT))*SPHP
      CALL DUDBRB(0,0,STHETA,0.D0,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +  DBCMS(3)/DBCMS(4))
c MOMENTA OF POMERON
      K(NIA1+1,1)=21
      K(NIA1+1,2)=KINT(2,2)
      K(NIA1+1,3)=2
      DO 10  KI=1,4
         P(NIA1+1,KI)=P(2,KI)-P(NPFIN,KI)
   10 CONTINUE
      P(NIA1+1,5)=-sqrt(ABS(DOT1(NIA1+1,NIA1+1)))
	Nia2 = Nia1 + 1
      T2GKI=-SNGL(P(NIA1+1,5)**2)
c      IF(IDEBUG.EQ.1) THEN
c         write(6,*) ' PARTDF: T2GKI,T2 ',T2GKI,T2
c         write(6,*) ' partdf: yx=',yx,dot1(nia1,2)/dot1(1,2)
c         write(6,*) ' partdf: xr=',xr,dot1(nia1,nia1+1)/dot1(nia1,2)
c      ENDIF
c for construction of parton in pomeron goto gamma pomeron system
c     write(6,*) ' back in gamma p system '
c     call dulist(1)

      DBCMS(1)=  P(NIA1,1) + P(NIA1+1,1)
      DBCMS(2)=  P(NIA1,2) + P(NIA1+1,2)
      DBCMS(3)=  P(NIA1,3) + P(NIA1+1,3)
      DBCMS(4)=  P(NIA1,4) + P(NIA1+1,4)
      spom = dot(dbcms,dbcms)
      shh = spom
c      write(6,*) ' rgsatrev spom,t2,pt ',spom,t2,pt
      IF(SPOM.LT.0.0D0) GOTO 30
      BOCHCK = (DBCMS(1)/DBCMS(4))**2 + (DBCMS(2)/DBCMS(4))**2
     +              + (DBCMS(3)/DBCMS(4))**2
      BOCHCK = DSQRT(BOCHCK)
      IF(BOCHCK.GT.0.99999999D0) goto 40
c      write(6,*) ' now goto gamma pomeron '
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  -DBCMS(1)/DBCMS(4),-DBCMS(2)/DBCMS(4),
     +  -DBCMS(3)/DBCMS(4))
      SPHI = DLANGL(P(NIA1,1),P(NIA1,2))
      CALL DUDBRB(0,0,0.D0,-sphi,0.d0,0.d0,0.d0)
      STHETA = DLANGL(P(NIA1,3),P(NIA1,1))
      CALL DUDBRB(0,0,-STHETA,0.d0,0.d0,0.d0,0.d0)
c end new


c  here construct the full event record in gamma pom system
c     call dulist(1)
c  reconstruct  four vectors of all partons
      Call PartKine
c     write(6,*) 'before rotate back'
c     call dulist(1)

      CALL DUDBRB(0,0,STHETA,0.d0,0.d0,0.d0,0.d0)
c     write(6,*) 'after theta rotate back'
c     call dulist(1)
      CALL DUDBRB(0,0,0.D0,sphi,0.d0,0.d0,0.d0)
c     write(6,*) 'after phi rotate back'
c     call dulist(1)
      CALL DUDBRB(0,0,0.D0,0.D0,
     +  DBCMS(1)/DBCMS(4),DBCMS(2)/DBCMS(4),
     +  DBCMS(3)/DBCMS(4))
c     write(6,*) 'after boost back'
c     call dulist(1)

      Return
   20 write(6,*) ' BOCHCK 1 ',BOCHCK
c     call dulist(1)
      Return
   30 write(6,*) ' SPOM ',spom
c     call dulist(1)
      Return
   40 write(6,*) ' BOCHCK 2 ',BOCHCK,ybj
c      call dulist(1)
      Return
      End
