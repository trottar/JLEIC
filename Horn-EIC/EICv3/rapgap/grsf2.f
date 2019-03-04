
      double precision function grsf2(x,s,alp,bet,a,b,ga,gb,gc,gd,
     #                                ge,gep)
      implicit double precision (a-z)
      grsf2=(s*x**a*(ga+gb*sqrt(x)+gc*x**b)+
     #      s**alp*exp(-ge+sqrt(gep*s**bet*log(1.d0/x))))*
     #      (1.d0-x)**gd
      return
      end
