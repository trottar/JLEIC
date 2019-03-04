

      function gamma(q)
      DOUBLE PRECISION   gamma,q
      DOUBLE PRECISION   aem,alphasl,pi,charge
      common   /cons/   aem,pi
comment on aem: since we use f2 aem is cancelled, so it doesnt matter
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      charge=1.D0/3.D0*1.D0/3.D0+2.D0/3.D0*2.D0/3.D0+1.D0/3.D0*1.D0/3.D0
      gamma=aem*alphasl(q)**2*pi**2*charge
      return
      end
