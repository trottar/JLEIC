*CMZ :  2.07/03 24/05/99  17.34.05  by  Hannes Jung
*-- Author :
      function dhformf(t)
	Implicit None
      DOUBLE PRECISION   dhformf
      DOUBLE PRECISION   t,dl
c

c use donnachie landshoff formfactor
c        f(t)=(4-2.8*t)/(4-t)*1/(1-t/0.7)^2
      dl =(4.d0-2.8d0*t)/(4.d0-t)* 1.d0/(1.d0-t/0.7d0)**2
c         write(6,*) ' dhformf : ',dl,t
      dhformf = dl**2
      return
      end
