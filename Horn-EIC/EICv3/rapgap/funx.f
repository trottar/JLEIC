*CMZ :  2.08/04 22/12/99  15.39.26  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.24.13  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  11.50.49  by  Hannes Jung
*-- Author :    Hannes Jung   22/01/95
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C---INTERFACE FOR CALLS OF PARTON DISTRIBUTIONS FROM HERACLES
C---TAKEN FROM PYSTFU (4.1.91) HS
C    changed for the use with pomeron structure function RASTFU
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FUNX(XINT)
	Implicit None
      DOUBLE PRECISION T2MIN,T2,XR,MP,XP2T,XX,WTDIST,XMAX
      REAL XPQ(-6:6),draprn
      DOUBLE PRECISION WTTEST,T2MN,T2MAX1,WTPO
      DOUBLE PRECISION XMAXV,XMINV,VMIN,VMAX,SMALL
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEND.
	Integer KPART,IVM
      REAL SCALE
      COMMON/PINT/ SCALE,KPART
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      COMMON/VMESON/IVM
	Real FUNX,WGMAX,XINT,ULMASS
	Integer I,KF
      EXTERNAL ULMASS
      DATA MP/0.938/
      DATA SMALL/1.d-3/
      XMAX=1.d0-XF
      FUNX = 0.
C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron
C ... XR  = E_pomeron/E_proton
ctest      XR = DBLE(XINT)
      XR = DBLE(XPR)*(XMAX/DBLE(XPR))**DBLE(XINT)
      IF(IVM.GE.1) THEN
         IF(IVM.GE.1.AND.IVM.LT.443) THEN
            VMAX = DBLE(1.020 + 2.*ULMASS(211))
            VMIN = DBLE(0.780 - 2.*ULMASS(211))
         ELSEIF(IVM.EQ.443) THEN
            VMAX = 2.D0*DBLE(ULMASS(421))
            VMIN = 2.D0*DBLE(ULMASS(4))
         ELSEIF(IVM.EQ.553) THEN
            VMAX = 2.D0*DBLE(ULMASS(521))
            VMIN = 2.D0*DBLE(ULMASS(5))
         ENDIF
         VMIN = VMIN+SMALL
         VMAX = VMAX-SMALL
c                  write(6,*) 'pomstr :',XMAX,XMIN
         XMAXV = DBLE(XPR)*(1 + VMAX**2/DBLE(SCALE))
         XMINV = DBLE(XPR)*(1 + VMIN**2/DBLE(SCALE))
c                  write(6,*) 'pomstr:',XMAXV,XMINV
         XR = XMINV*(XMAXV/XMINV)**DBLE(XINT)
         IF(XMAXV.GT.XMAX) RETURN
      ENDIF
      XP2T = DBLE(XPR)
      XX = XP2T/XR
      FUNX = 0.
      IF(XX.GE.0.99999) RETURN
      T2 = -99999.D0
      T2MIN=MP*MP*XR*XR/(1.d0-XR)
      T2MN = -(DBLE(SCALE) *(1.D0-XR/DBLE(XPR))
     +         +4.D0*DBLE(ULMASS(211))**2)
      IF(T2MN.LT.T2MAX) THEN
         T2MAX1 = T2MN
      ELSE
         T2MAX1 = T2MAX
      ENDIF

      T2 = T2MIN*((T2MAX1/T2MIN)**draprn())
c      write(6,*) ' FUNX: XR,XPR,XX',XR,XPR,XX
      IF(T2MIN.GE.T2MAX1) RETURN
      IF(XP2T.GE.XMAX) RETURN
      IF(XP2T.GE.XR) RETURN
c      write(6,*) ' pomstr ',NG,NPOM
      IF(NG.EQ.20.AND.NPOM.EQ.20) THEN
         KF=211
      ELSEIF(NG.EQ.21.AND.NPOM.EQ.21) THEN
         KF=111
      ELSE
         KF=100
      ENDIF
c      if(KF.EQ.211) THEN
c        write(6,*) ' FUNX: xx,XPR,XR ',XX,XPR,XR
c        write(6,*) ' FUNX: T2min = ',T2MIN,'T2MAX ',T2MAX
c        FUNX = 0
c        RETURN
c        endif
c      write(6,*) ' FUNX: XR,XP2T,T2mIN,T2',XR,XP2T,T2MIN,T2
c this is called to get maximum of differential pomeron distribution
      CALL RAT2DI(KF,XR,-T2,WTTEST)
      CALL RAT2IN(KF,XR,XP2T,T2MIN,-T2,WTDIST)
      IF(XX.GE.1.) THEN
         write(6,*) ' funx error ',xx,xr,xp2t
      endif
      CALL RASTFU(KF,SNGL(XX),SCALE,XPQ)
c      write(6,*) xpq(kpart),xpq(-kpart),kpart,SNGL(XX)
      FUNX = 0.
c       write(6,*) ' FUNX: WTDIST',WTDIST

      FUNX = XPQ(KPART)*SNGL(WTDIST)
      IF(IVM.EQ.0) THEN
         WTPO = T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAX/XP2T)*WTTEST
         FUNX = FUNX * SNGL(XR) * ALOG(SNGL(XMAX)/XPR)
      ELSE
         FUNX = FUNX * SNGL(XR * DLOG(XMAXV/XMINV))
         WTPO = T2*DLOG(T2MAX/T2MIN)*XR*DLOG(XMAXV/XMINV)*WTTEST
      ENDIF
      WGMAX=0.
      DO 10  I=-6,6
         IF(I.EQ.0) GOTO 10
         WGMAX = MAX(WGMAX,XPQ(I))
   10 CONTINUE
      WTPO=WTPO*DBLE(WGMAX)
      IF(WTPO.GT.POM) THEN
         POM = WTPO
ccc         write(6,*) ' WTPO,xr,t2 ',WTPO,XR,T2
      ENDIF
      RETURN

      END
