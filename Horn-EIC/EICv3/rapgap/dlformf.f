

      function dlformf(m,q,xb,i)
      DOUBLE PRECISION   dlformf,m,xb,q,xp
      DOUBLE PRECISION   a,b1,b2
      DOUBLE PRECISION   dexpin
      DOUBLE PRECISION   s,wmax,wmin
      DOUBLE PRECISION   ss,beta,xpp,t,qq,alphap,dl
      common   /parameter/ ss,qq,beta,xpp,t
      integer  i
c
      common   /hera/   s,wmax,wmin
c
      xp=xb*(q+m)/q
      a=dlog(1.d0/xp)
c
      if(i.eq.1) then
         dlformf=0.30119d0
      else if(i.eq.2) then
         b1=dexpin(0.35d0*a)
         b2=dexp(0.35d0*a)
         dlformf=-0.00191*b1*b2*a**3-0.085699*a+0.04559*dexp(2.d0*a)*
     +   dexpin(2.d0*a)+0.00546*a**2-0.052477*a*dexp(2*a)* dexpin(2.d0*
     +   a)-0.02274*a*b1*b2+0.0245331*a**2*b1*b2+ 0.380649-0.04558633*
     +   b1*b2
      elseif(i.eq.0) then
c use donnachie landshoff formfactor
c        f(t)=(1/x)^(alpha'*t)*(4-2.8*t)/(4-t)*1/(1-t/0.7)^2
c        alpha' = 0.25
         alphap = 0.25d0
         dl =(1.d0/xp)**(alphap*t)*(4.d0-2.8d0*t)/(4.d0-t)*
     +   1.d0/(1.d0-t/0.7d0)**2
c         write(6,*) ' dlformf : ',dl
         dlformf = dl**2
      end if
      return
      end
