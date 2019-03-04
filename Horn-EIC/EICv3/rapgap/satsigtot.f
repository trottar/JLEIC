*CMZ :  2.08/06 31/03/2000  11.20.32  by  Hannes Jung
*CMZ :  2.08/05 30/03/2000  15.33.05  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      Subroutine SATSIGTOT(Q2IN,XBJIN,SIGOUT,ESIGOUT,WEIOUT)
c
c Compute total diffractive sigma's (i.e. integrated over all masses)
c in microbarn at Q2IN,XBJIN by MC computation
c Finds also the maximal weight
      IMPLICIT NONE

      Double Precision Q2IN,XBJIN,SIGOUT,ESIGOUT,WEIOUT,W2IN
      Double Precision SIGXX, ESIGXX, DSigMax, DUM, weimc
	DATA DUM/0.d0/
      Double Precision DM, DW, DQ, EPS

      Double Precision Q2min, Q2max, W2min, W2max, Mxmin, Mxmax
      COMMON /GDLIMIT/ Q2min, Q2max, W2min, W2max, Mxmin, Mxmax

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision    DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT
      Double Precision    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT
      COMMON /CCOMPCROS/  DLQQT, DLQQL, DLQQGT, DLQQGL, DLTOT,
     &                    DCQQT, DCQQL, DCQQGT, DCQQGL, DLCTOT


      Integer         NQUAD
      COMMON /CNQUAD/ NQUAD

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGPARAM.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.
	Double Precision XMAX
      Integer NMaxWe, Iev
      Logical FIRST
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         FIRST=.FALSE.

         NQUAD = 24
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Parameters of Wuesthoff Model
         Call SATURPAR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ENDIF

      Stot = SSS
      SIGOUT = 0.
      ESIGOUT = 0.

      SIGXX = 0.
      ESIGXX = 0.
      DsigMax = 0.
      XMAX=1.d0-XF
	
      if(XBJIN.GE.XMAX) RETURN

      W2IN = Q2IN*(1.-XBJIN)/XBJIN

      Q2min = Q2IN
      Q2max = Q2IN+0.0005

      W2min = W2IN
      W2max = W2IN+0.0005
c
      Mxmin = 0.3
      Mxmax = max(Mxmin,dsqrt(Q2IN*(XMAX-XBJIN)/XBJIN))
c	write(6,*) ' Mx ',Mxmin,dsqrt(Q2IN*(XMAX-XBJIN)/XBJIN),Mxmax


      NMaxWe = 15000

      Do Iev =1,NMaxWe

         Call RANQ2W2Mx
         Call INPUTKINE

         Call COMPCROSS

         EPS = 2.*(1.-YBj)/(1.+(1.-YBj)**2)

         DLTOT = DLQQT+DLQQGT + EPS*(DLQQL)
         DLCTOT= DCQQT+DCQQGT + EPS*(DCQQL+DCQQGL) + DLTOT

         weimc = DLCTOT/Float(NMaxWe)

         SIGXX = SIGXX + weimc
         ESIGXX = ESIGXX + weimc**2

         if(DLCTOT.gt.DSigMax) Then
            DSigMax = DLCTOT
c            Write(6,*) 'sigtot ',DsigMax,ajac
         Endif

      Enddo
c
c  End of diffractive cross section computation by MC
c

C normalize cross sections

      SIGout = SIGXX/ 1000.
      ESIGout = sqrt(ESIGXX)/1000.

      WEIOUT = DsigMax

c      Write(*,10000) Q2IN,XBJIN,SIGOUT, ESIGOUT,WEIOUT,ajac
10000 format('Q2,X, SIGOUT, WEIOUT,ajac',6E15.4)

      Return
      End
