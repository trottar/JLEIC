      SUBROUTINE SIGMD(pt2,sigt,sigl)
c note that sigt,sigl is here defined as q2/4/pi*sigt (sigl)
c is actually the structure funtion
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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEND.
      DOUBLE PRECISION lambda2,mq,mp,mx2,beta0,mu0
      DOUBLE PRECISION PT2MD,PHIMD
      COMMON/PT/PT2MD,PHIMD
      COMMON/HARDPOM/PT2GEN,PHIGEN
      COMMON /PARMD/ beta,xpom,q2,t2,mq,eq2
      COMMON /SEL/IPHI
      REAL ULMASS
      EXTERNAL ULMASS,LUCHGE
      sigt = 0.d0
      sigl = 0.d0
      mp = 0.938d0
      mp = 1.D0
      mx2 = Q2*(1.d0-beta)/beta + t2
      alpha_p=1.08d0 + t2/2.d0
      f1t = (4.d0*mp**2 - 2.8d0*t2)/(4.d0*mp**2 - t2)
      f1t = f1t/(1.d0 - t2/0.7d0)**2
c      write(6,*) ' F2MD t2: ',t2
      beta0 = 2.d0
      mu0 = 1.1d0
      betamu = beta0*mu0
      alphas0 = 1.d0
      pi=2.d0*dasin(1.d0)
c      write(6,*) 'sigMD: beta,xpom,q2,t2,f1t',beta,xpom,q2,t2,f1t
c      write(6,*) ' pt2 ',pt2
      pl2 = mx2/4.d0 - mq**2 - pt2
      if(pl2.le.0.0D0.or.pt2.le.0.0d0) return
      pl = dsqrt(pl2)
      pt = dsqrt(pt2)

c now F2MD
      lambda2 = (pt**2 + mq**2)/(1.d0 - beta)
      alpha_s=alphas(dsqrt(lambda2))
      if(alpha_s.lt.0.d0) alpha_s=0.d0
      if(alpha_s.ge.1.d0) alpha_s=1.d0
      xm = mx2/(pt2+mq**2)
      xp = pt2/(pt2+mq**2)
      aa = (1.d0-beta)*xp
      sigpp = xm**2*(1.d0 -2.d0/xm)*xp*(1.d0-aa)**2
      sigpp = sigpp+ 0.25d0*xm**2*mq**2/(pt2+mq**2)*(1.d0-2.d0* aa)**2
      sigl = Q2/mx2 *xm *(1.d0 -2.d0*aa)**2
      sigpm = 2.d0*xm*xp*(1.d0-aa)**2
      sigp0 = -2.d0*dsqrt(q2)/dsqrt(mx2)*xm*pl/pt*xp
      sigp0 = sigp0*(1.d0-aa)*(1.d0-2.d0*aa)
      sigt = sigpp
      comfac = 1.d0/dsqrt(mx2)*f1t**2 *xpom**(2.d0* (1.d0-
     +alpha_p))/xpom*27.d0/2.d0/pi**3 *eq2*alpha_s/alphas0*betamu**4/
     +(mx2+q2)**2 *(1.d0-beta)
c            write(6,*) 'SIGMD: pt2,comfac,sigt',pt2,comfac,sigt
      epsilon =(1.d0 - dble(yy))/(1.d0 - dble(yy) + dble(yy)**2/2.d0)
      if(iphi.eq.1) then
         sigt = comfac/2.D0/pi*(sigt + epsilon*sigl
     +    - epsilon*dcos(2.d0*phimd)*sigpm
     +    - dsqrt(2.d0*epsilon*(1.d0+epsilon))*dcos(phimd)*sigp0)
         sigl = comfac*sigl
      else
         sigt = comfac*(sigt+epsilon*sigl)
         sigl = comfac*sigl
      endif
      if(sigt.lt.0.0) then
         write(6,*) ' sigt<0.0 ',comfac,epsilon,phimd,sigpp
         write(6,*) ' sigt1 ',- epsilon*dcos(2.d0*phigen)*sigpm
         write(6,*) ' sigt2 ',- dsqrt(2.d0*epsilon*(1.d0+epsilon))
         write(6,*) ' sigt3 ',dcos(phimd)*sigp0
      endif
      RETURN
      END
