*CMZ :  2.08/05 27/03/2000  16.12.39  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      SUBROUTINE RAN2BODY
*-- Author :   H. Kowalski
C   Routine generates the 2-body Kt variable and stores it in
C   the COMMON /C2BODVAR/  together with the corresponding Jacobian WT2
C
      Implicit NONE
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Integer QFKF
      Double Precision    MQUARK, MCHARM
      COMMON /CQMASS/     MQUARK, MCHARM, QFKF

      Integer QFKC, LUCOMP
      Double Precision    Kt2, WT2
      COMMON /C2BODVAR/   Kt2, WT2
      Double Precision   Kt2BL, Kt2BC
      COMMON /C2BOPT/    Kt2BL, Kt2BC

      Double Precision   Qpar2BL   ,QBpar2BL,   Qpar2BC   ,QBpar2BC
      Common /C2BOKINE/  Qpar2BL(4),QBpar2BL(4),Qpar2BC(4),QBpar2BC(4)

      Double Precision  Kt2Min, Kt2Max, Kt2JAC
      Double Precision  RANUNI, AALAM, MQU2, TWOPI, PHIRN, Kt,Kl2,Kl
      External          RANUNI, AALAM, LUCOMP


      Logical Debug2, Debug3


      MQU2  = MQUARK**2
c
c   2-parton phase space
c
c   create Kt2
c
      Kt2min = 0.0001
      Kt2max = AALAM(Mx2,MQU2,MQU2)/4./Mx2
c      Write(*,*) QFKF,MQUARK, MQU2, Mx2, Kt2Max

      Kt2 = DEXP(RANUNI(DLOG(Kt2min),DLOG(Kt2max)))
      Kt2JAC = Kt2*(DLOG(Kt2max/Kt2min))

      WT2 = Kt2JAC

      KL2 = Mx2/4.-Kt2-MQU2
      if(Kl2.lt.0.0) Then
         write(*,*) 'RAN2BODY- Kl2 lt 0', Kl2
         Kl2 = 0.
      Endif
      Kl = dsqrt(Kl2)
      Kt = dsqrt(Kt2)

      TWOPI = 2.*3.1415927
      PHIRN = RANUNI(0.D0,TWOPI)

      If(MQU2.lt.0.01) Then
         Qpar2BL(1) = Kt*DCOS(PHIRN)
         Qpar2BL(2) = Kt*DSIN(PHIRN)
         Qpar2BL(3) = -Kl
         Qpar2BL(4) = sqrt(Kt2+Kl2+MQU2)

         QBpar2BL(1) = -Kt*DCOS(PHIRN)
         QBpar2BL(2) = -Kt*DSIN(PHIRN)
         QBpar2BL(3) =  Kl
         QBpar2BL(4) = sqrt(Kt2+Kl2+MQU2)

         Kt2BL = dsqrt(Kt2)
      Endif

      If(MQU2.gt.1.) Then
         Qpar2BC(1) = Kt*DCOS(PHIRN)
         Qpar2BC(2) = Kt*DSIN(PHIRN)
         Qpar2BC(3) = -Kl
         Qpar2BC(4) = sqrt(Kt2+Kl2+MQU2)

         QBpar2BC(1) = -Kt*DCOS(PHIRN)
         QBpar2BC(2) = -Kt*DSIN(PHIRN)
         QBpar2BC(3) =  Kl
         QBpar2BC(4) = sqrt(Kt2+Kl2+MQU2)

         Kt2BC = dsqrt(Kt2)
      Endif


      Return
      End
