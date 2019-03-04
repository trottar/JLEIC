c---------------------------------------------------------------------
      subroutine xntotal(xnqt,xng,xnql,q2,x_pom,beta,v,z,t)
c
c
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION lim
c      double precision xnqt,fqt,xng,fg,xnql,fql,lim,
c     +                 q2,yb,beta,v,z,t,x_pom,jacob,pom,xpom,alp
c
c
      pi=2.0d0*asin(1.0d0)
      q02=1.0d0
      alp=0.25d0
      xng = 0.d0
      xnqt = 0.d0
      xnql = 0.d0
c
c
      xpom=pom(x_pom,-t,v)
      xpom = xpom/x_pom
      if (v.lt.q02/q2) then
         xng=0.0d0
      else
         xng=fg(z,v)/v
         xng=xng**2/q02*((1.0d0-beta/z)**2+(beta/z)**2)
         xng=xng*9.0d0/4.0d0*alp/8.0d0/pi*log(v*q2/q02)
c
         xng = xng/12.d0/z**2 * xpom * beta
c          write(6,*) ' xntotal xng,xpom',xng,xpom

      endif
c
      xnql=fql(beta,v)
      xnql=xnql**2/q2/v**3*4.0d0/3.0d0
c
      xnql = xnql/6.d0 * xpom * beta**3
c
c      write(6,*) ' xntotal xnql,xpom',xnql,xpom
      xnqt=fqt(beta,v)/v
      xnqt=xnqt**2/q02/6.0d0/(1.0d0-beta)
      xnqt=xnqt/sqrt(1.0d0-4.0d0*q02/q2*beta/v)
c add in new factor from ANL-HEP-PR-97-03
cccccccc      xnqt=xnqt*(1.0d0-2.0d0*q02/q2*beta/v)
      xnqt = xnqt/12.d0*beta *xpom
c      write(6,*) ' xntotal xnqt,xpom',xnqt,xpom
c
      lim=1.0d0-4.0d0*q02/q2*beta/v
      if (lim.lt.0.01d0) then
         xnql=xnql*20.0d0
         xnqt=xnqt*20.0d0
      else
         xnql=xnql/dsqrt(lim)
         xnqt=xnqt/dsqrt(lim)
      endif
c
      return
      end
