*CMZ :  2.07/03 24/05/99  16.51.32  by  Hannes Jung
*-- Author :    Hannes Jung   25/04/96
      SUBROUTINE USDIFFR(BETA,SCALE,XPQ,X_POM,T2)
	Implicit None
c note xpq is the parton density of the pomeron multiplied by the pomeron flux
c xpq = f_pomeron(x_pom,t) * beta * f(beta,scale)
      REAL T2,X_POM
      REAL BETA,SCALE
      REAL XPQ(-6:6)
c for example
      REAL beta02,rn2,epsp,alphp,pi,wtdist,xe,xp
	Integer I
      DATA beta02/58.74/,RN2/4./,EPSP/0.139/,alphp/0.3/,pi/3.1415/
c end example
      DO 10 I=-6,6
         XPQ(I)=0.0
   10 CONTINUE
c example : xq(x) = x(1-x) (only u,d quarks,gluons) and streng pomeron flux
C.... xq (X)= 6x(1-X)
c pomeron distribution from streng (hera proc. 1987)
      XE=RN2*ABS(T2)
      IF(XE.GT.170.0) THEN
         XE=170.
      ENDIF
      WTDIST = EXP(-XE)
      IF(WTDIST.LT.1.E-20) WTDIST=0.0
      WTDIST = WTDIST*(X_POM**(1.0-2.0*(1.0 + EPSP + ALPHP*T2)))
      WTDIST = WTDIST*beta02/16.0/PI
C.... beta q (beta)= 6/4 beta(1-beta)
      XP = 6.*(1.-beta)*beta
      XPQ(0) = XP*WTDIST
      XPQ(1) = XP/4.*WTDIST
      XPQ(2) = XPQ(1)
      XPQ(-1) = XPQ(1)
      XPQ(-2) = XPQ(2)
c      write(6,*) ' USDIFFR xp,wtdist,beta,x_pom,t2,scale:'
c     & ,xp,wtdist,beta,x_pom,t2,scale

      RETURN

      END
