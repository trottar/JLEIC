
      FUNCTION FMD(X)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
*KEEP,RGPARAS.
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

*KEND.
      DOUBLE PRECISION mq,mx2
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON /PARMD/ beta,xpom,q2,t2,mq,eq2
      REAL ULMASS
      EXTERNAL ULMASS,LUCHGE
      REAL X,FMD
      DATA SMALL/1.D-6/
      pi=2.d0*dasin(1.d0)
      pt2min=pt2cut(12)
      mx2 = Q2*(1.d0-beta)/beta + t2
      FLQ = 0.D0
      F2Q = 0.D0
      FMD = 0.
c      write(6,*) 'F2MD: beta,xpom,q2,t2,f1t',beta,x_pom,q2,t2,f1t
      wchi = 2.d0*pi
c try 1/pt**4
      pt2max = mx2/4.d0 - 4.d0*mq**2
c         write(6,*) pt2min,pt2max
      if(pt2min.gt.pt2max) return
      rn2 = dble(x)
      pt2 = pt2max*pt2min/(pt2max + rn2*(pt2min-pt2max))
      pl2 = mx2/4.d0 - mq**2 - pt2
      wpt2 = pt2**2*(pt2max-pt2min)/pt2max/pt2min
c multply wpt2 with jacobian pt -> pl : pt/pl
c and account for divergence at high pt --> pl --> 0
      if(rn2.gt.0.998d0) return
      wpt2=wpt2*dsqrt(pt2)/dsqrt(pl2)
      call sigmd(pt2,sigt,sigl)
c now from sigma to F2, FL
      F2Q = sigt*Q2/4.d0/pi * wpt2
      FLQ = sigl*Q2/4.d0/pi * wpt2
      IF(F2Q.LE.0.0) THEN
         write(6,*) 'FMD<0 ',sigt,sigl,wpt2
      ENDIF
c      write(6,*) 'F2Q ',F2Q
      pt2md=pt2
      FMD = SNGL(F2Q)
      RETURN
      END
