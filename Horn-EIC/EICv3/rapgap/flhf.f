
      FUNCTION FLHF(X1,X,Q2)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      REAL X1,FLHF,XS,sscal
      REAL X,Q2,LUCHGE,ULMASS,XX,Q2X
      COMMON /HEAVYF/XX,Q2X,KPHF
      REAL XPQ(-6:6)
      REAL SNGL
      EXTERNAL LUCHGE,ULMASS
      pi=2.d0*dasin(1.d0)
      flhf = 0.0
      f2 = 0.0D0
      xt = 0.0D0
      AM = DBLE(ULMASS(KPHF))
      SMIN = 4.D0*AM**2
      CHIMIN=DBLE(X)*(1.D0 + SMIN/DBLE(Q2))
      CHIMIN = DMAX1(CHIMIN,DBLE(X))
      XPMIN = DBLE(X)
      XPMAX = DBLE(X)/CHIMIN
c      write(6,*) ' flhf x,q2',x,q2
c       write(6,*) ' flhf xpmax,xpmin',xpmax,xpmin
c       write(6,*) ' flhf: CHIMIN,x,smin,q2',CHIMIN,x,smin,q2
      if(xpmax.lt.xpmin) then
c       write(6,*) ' flhf xpmax,xpmin',xpmax,xpmin
c       write(6,*) ' flhf: CHIMIN,x,smin,q2',CHIMIN,x,smin,q2
         return
      endif
      IF(XPMIN.EQ.0.0D0) THEN
         write(6,*) ' flhf XPMIN CHIMIN,X ',CHIMIN,X,AM
      ENDIF
      XP = XPMIN * (XPMAX/XPMIN)**DBLE(X1)
cc      XP = XPMIN + DBLE(X1)*(XPMAX-XPMIN)
      IF((4.d0*am**2*xp/DBLE(Q2)/(1.d0-xp)).GE.1.d0) RETURN
      beta2 = 1.d0 - 4.d0*am**2*xp/DBLE(Q2)/(1.d0-xp)
      beta = dsqrt(beta2)
      ZQMAX = 1.d0 + beta
      ZQMIN = 1.d0 - beta
      xksi = am**2/DBLE(q2)
c         IF(ZQMIN.GT.0.5) RETURN
c        write(6,*) ' flhf CHIMIN,X,am,kphf ',CHIMIN,X,AM,KPHF
c        write(6,*) ' flhf beta,zqmax,zqmin',beta,zqmax,zqmin,xp
      zlog = dlog(zqmax/zqmin)
      SCAL=4.D0*AM**2
      ALPHA_S=ALPHAS(DSQRT(SCAL))
c      write(6,*) ' flhf alphas ',alpha_s
      xg = dble(x)/xp
      xt = zlog*(0.5d0*xp - xp**2 *(1.d0-xp) + 2.d0*xksi*xp**2*(1.d0-
     +3.d0*xp) - 4.d0*xksi**2*xp**3) + beta*(4.d0*xp**2*(1.d0-xp) -
     +0.5d0*xp -2.d0*xp**2*xksi*(1.d0-xp))
c note the following line is different in massless case from Peccei/Ruckl
c    +        + beta*(3.d0*xp**2*(1.d0-xp) - 0.5d0*xp ! Peccei/Ruckl
      xl = 2.d0*beta*xp**2*(1.d0-xp) -4.d0*xp**3*xksi*zlog
      f2 = alpha_s * xt /pi
      fl = alpha_s * xl /pi
c      write(6,*) 'x,xp,xglu ',x,xp,xg
      f2 = f2/xp
      fl = fl/xp
      if(xg.ge.1.) write(6,*) ' flhf  xg>1 ',xg,x,xp,x1
      XS = SNGL(XG)
      sscal = sngl(scal)
      CALL RYSTFU(2212,XS,SSCAL,XPQ)
c      write(6,*) ' flhf xs,scal,xpq(0) ',xs,scal,xpq(0)
      wtg = dble(XPQ(0))
      f2 = f2 * wtg
      if(x1.gt.1..or.x1.lt.0.) THEN
         write(6,*) ' flhf: x1 = ',x1,AM
      endif
c      write(6,*) ' f2,xp ',f2,xp
c      write(6,*) 'flhf: x1 = ',x1,AM
c      write(6,*) ' flhf ',f2,xp
      FLHF = SNGL(fl* XP*DLOG(XPMAX/XPMIN))
cc      flhf = SNGL(f2* (XPMAX-XPMIN))
      RETURN
      END
