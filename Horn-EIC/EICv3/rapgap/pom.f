c-----------------------------------------------------
      double precision function pom(xp,t,v)
      implicit DOUBLE PRECISION (a-g,o-z)

c
c
      q02=1.0d0
      pm=1.09d0
      eps=0.08d0
      if (log(q02/v).gt.-2.0d0) then
         pm=pm/(log(q02/v)+3.0d0)**0.558d0
         eps=eps+0.0997d0*log(log(q02/v)+3.0d0)
      endif
      pm=pm**2
      eps=2.0d0*eps
      pm=pm*(0.05d0/xp)**eps
c
      fm=(4.0d0+2.8d0*t)**2/(4.0d0+t)**2
      form=fm/(1.0d0+t/0.7d0)**4*xp**(0.5d0*t)
c      form=exp(-2.0d0*b*t)*xp**(0.5d0*t)
c
      pom=pm*form
c      write(6,*) 'pom: xp,t,v ',xp,t,v
c      write(6,*) 'pom: pm,form',pm,form
c      write(6,*) ' old pom = ',pm*form
C new parameters from ANL-HEP-PR 97-03
c new pomeron flux (gluon density)
      A = 0.877d0
      B = 0.133d0
      C = 0.596d0
      x0 = 0.05d0
      alpom = 0.085
      scale = 1.d0/v
      lnscale=log(scale/q02)
      if(lnscale.lt.0d0) lnscale=0.d0
      if(scale.gt.q02) then
         alpom = alpom + B*log(lnscale+1.d0)
      endif
c here we use (x0/xp)**alpom instead of (x0/xp)**(1-alpom) as in the paper
c which is presumably a typing error
      pm = A * (x0/xp)**(alpom) *(lnscale + 1.d0)**(-C)
c new form of t dependence
      B_D = 5.9d0
      form = exp(-B_D*t)
c      form = form/B_D
c      write(6,*) ' wuesthoff v,t,xp:',v,t,xp,pm,form
      pom = pm**2 * form
c      write(6,*) ' new pom = ',pm*form

c
      return
      end
