c---------------------------------------------------------------
      double precision function fg(z,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION ng,lnln

c
c      double precision z,v,root,lnln,ng,xa,xb
c
      if ((1.0d0-z)*v**2/(v-1.0d0)**2.lt.1.0d-3
     +                          .and.v.gt.1.0d0) then
         ng=(3.0d0*v-5.0d0)*v/(v-1.0d0)**2
         ng=ng+2.0d0*(1.0d0+1.0d0/(v-1.0d0)**3)*log(v)
         ng=ng*(1.0d0-z)**2
      elseif (v.gt.1.0d3) then
         ng=2.0d0*(z-1.0d0-z**2*log(z)-(1.0d0-z)**2*log(v))
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*z)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*z+v)
         lnln=lnln/(root+(1.0d0-2.0d0*z)*v+1.0d0)
         lnln=log(lnln)
         xa=1.0d0+v-2.0d0*z*(1.0d0-z)
         xb=root-2.0d0*z*(1.0d0-z)/root
c
         ng=xa*log(1.0d0/z)
         ng=ng+(2.0d0*v-xa+xb)*log(v)
         ng=ng+xb*lnln
      endif
      ng=ng/(1.0d0-z)
c       write(6,*)'ng= ',ng,' v = ',v,' beta = ',beta
c
      fg=ng
      return
      end
