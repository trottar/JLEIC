c---------------------------------------------------------------
      double precision function fqt(beta,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION nqt,lnln
c
c      double precision beta,v,root,lnln,nqt,xa,xb
c
      if (v.gt.1.0d3) then
         nqt=2.0d0*((1.0d0-beta)*log(v)-beta*log(beta))
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*beta)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*beta+v)
         lnln=lnln/(root+(1.0d0-2.0d0*beta)*v+1.0d0)
         lnln=log(lnln)
         xa=1.0d0-2.0d0*beta
         xb=(xa+v)/root
c
         nqt=xa*log(1.0d0/beta)
         nqt=nqt+xb*lnln
         nqt=nqt+(xb-xa)*log(v)
      endif
c
      fqt=nqt
c      write(6,*)'nqt=',nqt,' v = ',v,' beta = ',beta
      return
      end
