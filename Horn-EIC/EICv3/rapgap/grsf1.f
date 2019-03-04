
c
      double precision function grsf1(x,s,alp,bet,a,b,ga,gb,gc,gd,
     #                                ge,gep)
      implicit double precision (a-z)
      grsf1=(x**a*(ga+gb*sqrt(x)+gc*x**b)+
     #      s**alp*exp(-ge+sqrt(gep*s**bet*log(1.d0/x))))*
     #      (1.d0-x)**gd
      return
      end
