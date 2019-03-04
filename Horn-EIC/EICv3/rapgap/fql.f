c------------------------------------------------------------
      double precision function fql(beta,v)
      implicit DOUBLE PRECISION (a-g,o-z)
      DOUBLE PRECISION nql,lnln

c
c      double precision beta,v,root,lnln,nql
c
      if (v.gt.1.0d3) then
         nql=log(beta*v)
      else
         root=sqrt(v**2+2.0d0*(1.0d0-2.0d0*beta)*v+1.0d0)
c
         lnln=root-(1.0d0-2.0d0*beta+v)
         lnln=lnln/(root+(1.0d0-2.0d0*beta)*v+1.0d0)
         lnln=log(lnln)
c
         nql=log(1.0d0/beta)
         nql=nql+lnln/root
         nql=nql+(1.0d0/root-1.0d0)*log(v)
      endif
c
      fql=nql
c       write(6,*)'nql=',nql,'v = ',v,' beta = ',beta
      return
      end
