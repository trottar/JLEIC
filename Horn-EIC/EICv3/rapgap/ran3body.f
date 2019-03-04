*CMZ :  2.08/05 27/03/2000  16.13.23  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      SUBROUTINE RAN3BODY
*-- Author :   H. Kowalski
C   Routine generates the 3-body variables KtG2, QtG2, SM2, ZW, KW2
C   and stores it in the COMMON /C3BODVAR/  together with
C   the corresponding Jacobian WT3
C   The variable YCUT in COMMON /C3BODVAR/ is set in the main program
C   it is necessary when the process QGQb is evaluated (not used in the
C   present version)
C
      Implicit NONE
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf
      Double Precision    PGAMCM   , PPOMCM,   PGPOMCM
      Common /CINPKINE/   PGAMCM(5), PPOMCM(5),PGPOMCM(5)

      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Double Precision    Kpar3B(4), Qpar3B(4) ,QBpar3B(4)

      Double Precision    Kpar3BL   , Qpar3BL    ,QBpar3BL
      Double Precision    Kpar3BC   , Qpar3BC    ,QBpar3BC
      Common /C3BOKINE/   Kpar3BL(4), Qpar3BL(4) ,QBpar3BL(4),
     &                    Kpar3BC(4), Qpar3BC(4) ,QBpar3BC(4)

      Double Precision    Kt2, Qt2, SM2, ZW, KW2, WT3, YCUT
      COMMON /C3BODVAR/   Kt2, Qt2, SM2, ZW, KW2, WT3, YCUT
      Double Precision    Kt3BL, Qt3BL, ZW3BL, Kt3BC, Qt3BC, ZW3BC
      COMMON /C3BOPT/     Kt3BL, Qt3BL, ZW3BL, Kt3BC, Qt3BC, ZW3BC

      Double Precision       AlphaK, BetaK, AlphaQ, BetaQ

      Double Precision  Kt2min,  Kt2max,  KtJAC
      Double Precision  ZWMin,   ZWMax,   ZWJAC, Mx
      Double Precision  oZWMin,   oZWMax
      Double Precision  Qt2Min, Qt2Max, QtJAC, Qt2dum
      Double Precision  Pt2max
      Double Precision  RANUNI, AALAM, MQU2, MQU
      External          RANUNI, AALAM

      Double Precision    PGAMpri(4), PKtQt(2), PQQB(4), MQQB2
      Double Precision    PPOMpri(4)
      Double Precision  TWOPI, PHIRN, PHIKtQt, XSIGN,SIGN
      Double Precision  Kt, Qt, KtQt2, SS, Spri
      Double Precision  Left, Right, Leftpri, FourDot
      Double Precision  AA, BB, CC, Qt2m, SQ

      Integer I

      Logical Debug

      DEBUG = .false.
c
c   3-parton phase space
c
c   create Kt2
c
      TWOPI = 2.*3.1415927
      Mx = dsqrt(Mx2)
      WT3 = 0.
      IF(Mx.lt.1.) Return

      Kt2min =  0.001
c      Kt2min =  0.1
cH.K. 19.12.99 introduce charm mass
c  define MQU2 properly
      MQU2  = MQUARK**2
      Kt2max = AALAM(Mx2,0.0,4.D0*MQU2)/4./Mx2
cH.K. 19.12.99
c      if(Kt2max.gt.Q2) Kt2max = Q2

      Kt2   = DEXP(RANUNI(DLOG(Kt2min),DLOG(Kt2max)))
      KtJAC = Kt2*(DLOG(Kt2max/Kt2min))
      if(KtJAC.lt.0.0) Then
          write(*,*) ' Ktjac lt 0', Mx2, MQU2,Kt2max
          KtJAC = 0.
      Endif
c
c   create Spri, ZW
c
      Spri = (Mx2+Tdf+Q2)/2.
c      write(*,*) 'Spri 1 ',Spri, Q2,Mx2,Kt2,Tdf
      Spri = Spri*(1.+dsqrt(1.-4.*Tdf*Q2/(Mx2+Tdf+Q2)**2))
c      write(*,*) 'Spri 2 ',Spri

cH.K. 19.12.99 introduce charm mass
c old      ZWmin = ((Spri-Q2)/2./Spri)*(1.-dsqrt(1.-4.*Kt2/Mx2))+Q2/Spri
c old      ZWmax = ((Spri-Q2)/2./Spri)*(1.+dsqrt(1.-4.*Kt2/Mx2))+Q2/Spri

      ZWmin = ((MX2+4.*MQU2)/2./(Spri-Tdf))
     &   *(1.-dsqrt(1.-4.*Mx2*(Kt2+4.*MQU2)/(Mx2+4.*MQU2)**2))+Q2/Spri
      ZWmax = ((MX2+4.*MQU2)/2./(Spri-Tdf))
     &   *(1.+dsqrt(1.-4.*Mx2*(Kt2+4.*MQU2)/(Mx2+4.*MQU2)**2))+Q2/Spri
cH.K. 19.12.99

      if(ZWmin.gt.ZWmax) then
        ZWMIN = ZWMAX
        write(*,*) 'ZWMIN ',ZWmin,ZWmax,Q2,Mx2,Kt2
        write(*,*) 'Spri  ', Spri,Tdf
      Endif
      if(ZWmin.lt.Beta) then
        write(*,*) 'ZWMIN lt Beta ',Q2,Mx2,Kt2
      Endif
      if(ZWmin.gt.1.0) then
        write(*,*) 'ZWMIN gt 1. '
      Endif

      ZW = DEXP(RANUNI(DLOG(ZWmin),DLOG(ZWmax)))
      ZWJAC = ZW*(DLOG(ZWmax/ZWmin))
      if(ZWJAC.lt.0.0) ZWJAC = 0.

      If((ZW*Spri-Q2).lt.0.) Then
        write(*,*) ' ZW*Spri-Q2 < 0 '
        WT3 = 0.
        Return
      EndIf
      SM2 = ZW*Spri-Q2

      If(SM2.gt.Mx2) Then
        if((Sm2-0.5).gt.Mx2) Then
          write(*,*) ' SM2 gt Mx2',SM2,Mx2,Q2,W2,Tdf
        Endif
        WT3 = 0.
        Return
      Endif


c
c   create BetaK, AlphaK
c
      BetaK = 1.-ZW
      AlphaK = Kt2/(BetaK*Spri)

      PhiKtQt = RANUNI(0.D0,TWOPI)

cH.K 19.12.99
c old       Qt2min = Kt2

c old      Qt2max = (dsqrt(SM2)/2.)*(1.-(Tdf/Spri)-AlphaK-Kt2/SM2)/
c old     &      (dsqrt(1.-(Tdf/Spri)-AlphaK)+dsqrt(Kt2/SM2)*DCOS(PhiKtQt))

c      Fac1 =  1.-(Tdf/Spri)-AlphaK-Kt2/SM2
c      Fac2 =  1.-(Tdf/Spri)-AlphaK-Kt2/SM2*(DCOS(PhiKtQt)**2)
c      Fac4 = dsqrt((1.-(Tdf/Spri)-AlphaK)*(1.-4.*MQU2/SM2*Fac2/Fac1**2))
c      Fac4 = Fac4-dsqrt(Kt2/SM2)*DCOS(PhiKtQt)
c      Qt2max = dsqrt(SM2)/2.*Fac1/Fac2*Fac4
c      Qt2max = Qt2max**2
cH.K 19.12.99
cH.K.17.01.00

      AA = 1.-AlphaK-(Tdf/Spri)
      BB = AA-(Kt2/SM2)
      CC = AA-(Kt2/SM2)*(DCOS(PhiKtQt)**2)

      SQ = dsqrt(AA*(1.-4.*MQU2/SM2*CC/BB**2))
      Qt2max= 0.5*dsqrt(SM2)*BB/CC*(SQ-dsqrt(Kt2/SM2)*DCOS(PhiKtQt))

      Qt2min=-0.5*dsqrt(SM2)*BB/CC*(SQ+dsqrt(Kt2/SM2)*DCOS(PhiKtQt))
      if(Qt2min.lt.0.0) Qt2min = 0.0

      Qt2max = Qt2max**2
      Qt2min = Qt2min**2

      Qt2m = Kt2-MQU2
      if(Qt2m.lt.0.001)  Qt2m  = 0.001

      if(Qt2min.lt.Qt2m) Qt2min= Qt2m

c      Write(*,1001)  Qt2min,Qt2max,ZW,BETA,Kt2,Spri,SM2,PhiKtQt
c      Write(*,1001)  SQ,AA,BB,CC, MQU2

      if(Qt2max.lt.Qt2min) Then
c this case may happens when the  mass of the qq system (SM2) is small
c in this case Qt cannot be larger then (roughly) SM2/4.
c         Write(*,1001)  Qt2min,Qt2max,ZW,BETA,Kt2,MQU2
       Qt2dum = (dsqrt(SM2)/2.)*(1.-(Tdf/Spri)-AlphaK-Kt2/SM2)/
     &      (dsqrt(1.-(Tdf/Spri)-AlphaK)+dsqrt(Kt2/SM2)*DCOS(PhiKtQt))
       Qt2dum = Qt2dum**2
c       Write(*,1001) Qt2m,Qt2dum,Mx2,Sm2, tdf

c      oZWmin = ((Spri-Q2)/2./Spri)*(1.-dsqrt(1.-4.*Kt2/Mx2))+Q2/Spri
c      oZWmax = ((Spri-Q2)/2./Spri)*(1.+dsqrt(1.-4.*Kt2/Mx2))+Q2/Spri
c      write(*,*) ozwmin, zwmin
c      write(*,*) ozwmax, zwmax


         WT3 = 0.
         Return
 1001   Format(' Qt2min  Qt2max',10F9.4)
      Endif

      if(Qt2min.le.1.D-8) Then
c        Write(*,*) ' Qt2min', Qt2min,ZW,BETA,Kt2
        Qt2min = 1.D-8
      Endif

      Qt2 = DEXP(RANUNI(DLOG(Qt2min),DLOG(Qt2max)))
      QtJAC = Qt2*(DLOG(Qt2max/Qt2min))
      if(QtJAC.lt.0.0) QtJAC = 0.

      KW2 = Kt2/(1.-ZW)

      PHIRN = RANUNI(0.D0,TWOPI)

      Kt = dsqrt(Kt2)
      Kpar3B(1) = Kt*DCOS(PHIRN)
      Kpar3B(2) = Kt*DSIN(PHIRN)

      PhiRN = PhiKtQt+PhiRN
      If(PhiRn.gt.TWOPI) PhiRN = PhiRN-TWOPI

      Qt = dsqrt(Qt2)
      Qpar3B(1) = Qt*DCOS(PHIRN)
      Qpar3B(2) = Qt*DSIN(PHIRN)

      PKtQt(1) = Kpar3B(1) + Qpar3B(1)
      PKtQt(2) = Kpar3B(2) + Qpar3B(2)

      KtQt2 = PKtQt(1)**2 + PKtQt(2)**2

      XSIGN = RANUNI(0.D0,1.D0)

                         SIGN =  1.D0
      If(XSIGN.le.0.5D0) SIGN = -1.D0

      SS = (1.-(Tdf/Spri)-AlphaK+(Qt2-KtQt2)/SM2)**2
     *   - 4.*(1.-(Tdf/Spri)-AlphaK)*(Qt2+MQU2)/SM2
      If(SS.lt.0.) Then
c      write(*,*) ' '
c      write(*,*) '  SS lt 0', SS, MQU2, Kt2,Qt2, SM2
      WT3 = 0.
      return
      Endif

      AlphaQ =0.5*(1.-(Tdf/Spri)-AlphaK+(Qt2-KtQt2)/SM2+SIGN*dsqrt(SS))
      BetaQ  = Qt2/AlphaQ/Spri

      do i=1,4
       PGAMpri(i)=(PGAMCM(i)+  Q2*PPOMCM(i)/Spri)/(1.-Tdf*Q2/(Spri)**2)
       PPOMpri(i)=(PPOMCM(i)+ Tdf*PGAMCM(i)/Spri)/(1.-Tdf*Q2/(Spri)**2)
      EndDo

      Do i=3,4
        Kpar3B(i) = AlphaK*PGAMpri(i) + BetaK*PPOMpri(i)
        Qpar3B(i) = AlphaQ*PGAMpri(i) + BetaQ*PPOMpri(i)
      EndDo
      if((Qpar3B(4).lt.0.).or.(Kpar3B(4).lt.0.)) Then
       write(*,*) 'Energies lt 0'
      Endif

      do i=1,4
       QBpar3B(i) = PGAMCM(i)+PPOMCM(i)-Qpar3B(i)-Kpar3B(i)
       PQQB(i)    = QBpar3B(i) + Qpar3B(i)
      EndDo

      MQQB2 = FourDot(PQQB,PQQB)

c      write(*,*) TDF, FourDot(QBpar3B,QBpar3B),FourDot(Qpar3B,Qpar3B),
c     *           FourDot(Kpar3B,Kpar3B)

      If(SM2.lt.MQQB2) Then
       Write(*,*) ' Sm2 too small ', MQQB2, SM2
       SM2 = 0.
      Endif

      WT3 = KtJAC*ZWJAC*QtJAC/dsqrt(SS)

      If(MQU2.lt.0.01) Then
         Do i=1,4
            Kpar3BL(i)  = Kpar3B(i)
            Qpar3BL(i)  = Qpar3B(i)
            QBpar3BL(i) = QBpar3B(i)
         EndDo

         Kt3BL = dsqrt(Kt2)
         Qt3BL = dsqrt(Qt2)
         ZW3BL = ZW
      Endif

      If(MQU2.gt.1.00) Then
         Do i=1,4
            Kpar3BC(i)  = Kpar3B(i)
            Qpar3BC(i)  = Qpar3B(i)
            QBpar3BC(i) = QBpar3B(i)
         EndDo

         Kt3BC = dsqrt(Kt2)
         Qt3BC = dsqrt(Qt2)
         ZW3BC = ZW
      Endif

c      if(DEBUG) Then
c      write(*,*) ' '
c      write(*,*) '    SS,  AlphaK,   Kt2,    Qt2,    KtQt2,   SM2,  ZW'
c       write(*,1002)  SS, AlphaK,Kt2,Qt2,KtQt2,SM2,ZW,BETA
 1002  Format(6f8.2,2f8.3)
c      write(*,*) ' '
c      write(*,*) 'Wt3',ZWJAC,QtJAC,KtJAC
c      Endif

      Return
      End
