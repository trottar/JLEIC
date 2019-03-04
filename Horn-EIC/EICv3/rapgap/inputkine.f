*CMZ :  2.08/05 27/03/2000  16.11.42  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      Subroutine INPUTKINE
      IMPLICIT NONE
c
c   author    H.Kowalski
c   reconstruction of input parton (gamma, pommeron) four vectors from
c   dynamic (generated) variables  Q2,Mx2,W2
c   results is stored in Common /CINPKINE/
c
      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf
      Double Precision    PGAMCM   , PPOMCM,   PGPOMCM
      Common /CINPKINE/   PGAMCM(5), PPOMCM(5),PGPOMCM(5)

      Double Precision   MQU2, MQU
      DATA MQU /0.01/, MQU2 /0.0001/
      Double Precision  DUM, FourDot, RANUNI
      Double Precision   TWOPI
      Double Precision   Egam, Pgam, Epom, Mx
      EXTERNAL RANUNI

      Integer I

      LOGICAL DEBUG

      DEBUG = .false.
      TWOPI = 2.*3.1415927
c
c     Computation of gamma* and pomeron 4-vectors
c     in gamma-pomeron CMS
c
      Mx = dsqrt(Mx2)
      Egam = (Mx2-Q2+Tdf)/2./Mx
      Epom = (Mx2+Q2-Tdf)/2./Mx
      Pgam = ((Mx2+Q2+Tdf)**2-4.*Q2*Tdf)/4./Mx2
      If(Pgam.lt. 0.) write(*,*) '**** WARNING Pgam < 0  ****'
      Pgam = Dsqrt(Pgam)


      PGAMCM(1) = 0.
      PGAMCM(2) = 0.
      PGAMCM(3) = -Pgam
      PGAMCM(4) =  Egam

      PPOMCM(1) = 0.
      PPOMCM(2) = 0.
      PPOMCM(3) = Pgam
      PPOMCM(4) = Epom

      PGAMCM(5) = -Q2
      PPOMCM(5) = -Tdf

      Do i=1,4
         PGPOMCM(i)=PGAMCM(i)+PPOMCM(i)
      Enddo
      PGPOMCM(5) = dsqrt(FourDot(PGPOMCM,PGPOMCM))

      If(DEBUG) THEN
         Write (*,10000) Q2, sqrt(W2), sqrt(Mx2), Tdf
10000 Format('INPKINE- Q2, W, Mx Tdf',4F8.1/)

         Write (*,10100) (PGAMCM(i), i=1,5)
         Write (*,10100) (PPOMCM(i), i=1,5)
         Write (*,10100) (PGPOMCM(i), i=1,4)
10100 Format('   ',5F10.3)
         DUM = FourDot(PGPOMCM,PGPOMCM)
         Write(*,*) 'Mx  ', sqrt(dum)
         Write(*,*) ' '
      EndIf

      Return
      End
