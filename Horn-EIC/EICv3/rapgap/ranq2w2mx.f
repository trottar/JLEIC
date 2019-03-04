*CMZ :  2.08/05 27/03/2000  16.11.13  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      SUBROUTINE RANQ2W2MX
*-- Author :   H. Kowalski
C     Routine generates variables (Q2,W2,Mx2) and the
C     corresponding Jacobian factor AJAC in Common /CQ2W2MX/
C     The variables (Q2,W2,Mx2) are generated between limits
c     given in COMMON /GDLIMIT/
C     The routine computes also from Q2,W2,Mx2 the  variables
c     YBj, XBj, Xpom and Beta and puts them in COMMON /GDVARB1/
c
c     it generates also the variable Tdf (which plays a role of t)
c     Tdf generation is preliminary because it is to primitive
C-----------------------------------------------------------------------
      Implicit NONE

      Double Precision Q2min, Q2max, W2min, W2max, Mxmin, Mxmax
      COMMON /GDLIMIT/ Q2min, Q2max, W2min, W2max, Mxmin, Mxmax

      Double Precision  Stot, Q2, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2, W2, Mx2, AJAC
      Double Precision  YBj, XBj, Xpom, Beta, Tdf
      COMMON /GDVARB1/  YBj, XBj, Xpom, Beta, Tdf

      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

      DOUBLE PRECISION Q2JAC, Yjac, MxJAC

      Double Precision  RANUNI
      External          RANUNI
      Double Precision  W2min1, W2max1,Mx2min,Mx2max
      Double Precision  W, Mx, Rsta, Rend, Tmin, Tend, mprot

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.


C
C- generate the Q2 variable
C-
        Q2 = 	Q2MIN
	Q2JAC = 1.	
C
C- generate the W2 and YBj variable
C-
        W2min1 = W2min
        If(W2MIN1.lt.Q2) W2MIN1 = Q2
C-
        W2 =   W2MIN1
	YJAC = 1.
C-
        Ybj = (W2+Q2)/Stot
C-
C- generate the Mx variable
C-

      Mx = DEXP(RANUNI(DLOG(Mxmin),DLOG(Mxmax)))
      MxJAC = Mx * DLOG(Mxmax/Mxmin)

      Mx2 = Mx**2

      AJAC = Q2jac*Yjac*Mxjac

c  Generate the (positive) t variable in the range t=Tmin to t= Tend
c      Tend = 5.
      Tend = T2MAX
      If(Tend.Gt.Q2) Tend = Q2 - 0.1
      Tend = -Tend
      Tmin = 0.
      Rsta = DEXP(Tend)
      Rend = DEXP(Tmin)
      Tdf = -DLOG(RANUNI(Rsta,Rend))/6.

c       Tdf = 0.
C-
C-  compute the variables xb, xg, eta
C-
      Xbj = Q2/STOT/YBj
      beta = Q2/(Mx2 + Q2 + Tdf)
      Xpom = (Mx2 + Q2 + Tdf)/(W2 + Q2)

      R02 = (XPOM/X0)**LAM


c      W =   sqrt(W2)
c      write(*,1001) Q2,Mx,W,R02,Tdf

c  Tmin condition
      mprot = 0.9383
      tmin = mprot**2*Xpom**2/(1.-Xpom)

c      write(6,*) ' ajac =',ajac,Q2jac,Yjac,Mxjac,Mxmax,Mxmin
	If(Tdf.lt.tmin) Then
         AJAC = 0.
c         write(*,*) '***  t < Tmin  ***', Tdf, tmin
      Endif
C-


   10 RETURN
      END
