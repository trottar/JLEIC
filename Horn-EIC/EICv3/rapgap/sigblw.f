
      subroutine sigblw(pt,sigt,sigl)
c
c     <sig> berechnet den wq als funktion von pt^2, die anderen
c     kinematischen variablen werden im common uebergeben.
c     pt,q und m stehen immer fuer physikalisch pt^2,Q^2,M^2 .
c
c     die gluonstrukturfunktion wird im up grv_nl berechnet,
c     die ableitung der strfkt. in dgrv_nl.
c
c     das up gamma berechnet einen vorfaktor. es ruft ein
c     programm zur berechnung von running alpha_s auf.
c     durch das up dlformf wird implizit die t-integration
c     ausgefuehrt.
c
c     sigt, sigl, und sigi repraesentieren den gamma^* - p
c     wq fuer transversale und longitudinale photonen sowie
c     den interferenzterm.
c
      implicit none
      DOUBLE PRECISION   sigt,sigl,sigi,sigtt,epsilon
      DOUBLE PRECISION   dlformf
      DOUBLE PRECISION   s,q,xb,pt,m,beta,xp,t
      DOUBLE PRECISION   fact,facl,pi
      DOUBLE PRECISION   grv_nl,dgrv_nl
      DOUBLE PRECISION   gamma
      DOUBLE PRECISION   phi
c
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
      INTEGER IPHI
      COMMON /SEL/IPHI
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEND.
c
      pi=2.d0*dasin(1.d0)
      fact=4.d0/3.d0
      facl=4.d0/3.d0
c
      m=q/beta-q
      xb=beta*xp
c
c
      sigtt=gamma(pt*(q+m)/m)*
     +      fact*
     +     (1.d0-2.d0/m*pt)/dsqrt(1.d0-4.d0*pt/m)/pt**2*
     +     m**2/(q+m)**4*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*beta    +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       1.5d0*(1.d0-4.d0/3.d0*beta)  -
     +       1.d0*(1.d0-beta-(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )**2
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dx_pom
      sigtt = sigtt*q/xb
c
c
      sigl=gamma(pt*(q+m)/m)*
     +     facl*1.d0/dsqrt(1.d0-4.d0*pt/m)*
     +     q/(q+m)**4/pt*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*(2.d0*beta-1.d0) +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       3.d0*(1.d0-beta)  -
     +       (1.d0-beta-2.d0*(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )**2
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dbeta
      sigl = sigl*q/xb
c
      sigi=gamma(pt*(q+m)/m)
     +      /m**2*dsqrt(pt/q)*4.d0/3.d0*
     +      m*q*m**2/pt**2/(q+m)**4*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*beta    +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       1.5d0*(1.d0-4.d0/3.d0*beta)  -
     +       1.d0*(1.d0-beta-(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )*
     +     (
     +      grv_nl( xb*(q+m)/q,pt*(q+m)/m)*(2.d0*beta-1.d0) +
     +      1.d0*
     +      pt*dgrv_nl( xb*(q+m)/q,pt*(q+m)/m)*(
     +       3.d0*(1.d0-beta)  -
     +       (1.d0-beta-2.d0*(1.d0-beta)*dlog(beta)+dlog(beta))  )
     +     )
     +     *dlformf(m,q,xb,0)
c include jacobian dM**2 --> dbeta
      sigi = sigi*q/xb
      epsilon =(1.d0 - dble(yy))/(1.d0 - dble(yy) + dble(yy)**2/2.d0)
      if(iphi.eq.1) then
         sigt =    ((1.d0-2.d0*epsilon*dcos(2.d0*phi)*1.d0*
     +             pt/m/(1.d0-2.d0*pt/m)
     +              )*sigtt
     +             + (2.d0-dble(yy))*dsqrt(1.d0-dble(yy))*sigi*dcos(phi)
     +               /(1.d0 - dble(yy) + dble(yy)**2/2.d0)
     +             + epsilon*sigl
     +              )
     +              /2.d0/pi
         sigl = sigl/2.d0/pi
      else
         sigt = sigtt + epsilon*sigl
      endif
      if(sigt.lt.0.0) then
         write(6,*) ' sigt<0.0',sigtt,epsilon
         write(6,*) ' sigt1',-2.d0*epsilon*dcos(2.d0*phi)*1.d0*
     +    pt/m/(1.d0-2.d0*pt/m)
         write(6,*) ' sigt2',(2.d0-dble(yy))*dsqrt(1.d0-dble(yy))
     +    *sigi*dcos(phi)
         write(6,*) ' sigt3 yy',yy,q/xb/s,xb,phi
         write(6,*) ' sigl sigt',sigl,sigt
      endif
c
      return
      end
