      FUNCTION FXT2(XINT,T2INT)
	Implicit None
      DOUBLE PRECISION T2MIN,T2,XR,MP,XP2T,XX,WTDIST,XMAX,WTPI
      REAL XPQ(-6:6),WGMAX
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
	Real SCALE,ULMASS
	Integer KPART,IVM
      COMMON/PINT/ SCALE,KPART
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      COMMON/VMESON/IVM
	Integer I,KF
	Double Precision T2MAX1,T2MN
	Real FXT2,XINT,T2INT
      EXTERNAL ULMASS
      DATA MP/0.938/
      DATA SMALL/1.d-3/
      XMAX=1.d0-XF
C ... XP2 = E_glue /E_proton
C ... XX  = E_glue/E_pomeron
C ... XR  = E_pomeron/E_proton
c      XR = DBLE(XINT)
      FXT2 = 1.
c      write(6,*) NG,NPOM
      IF(NG.EQ.20.AND.NPOM.EQ.20) THEN
         KF=211
      ELSEIF(NG.EQ.21.AND.NPOM.EQ.21) THEN
         KF=111
      ELSE
         KF=100
      ENDIF
      XR = DBLE(XPR)*(XMAX/DBLE(XPR))**DBLE(XINT)
      FXT2 = SNGL(XR * DLOG(XMAX/DBLE(XPR)))
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
c                  write(6,*) 'pomstr:',XMAX,XMIN
         XMAXV = DBLE(XPR)*(1 + VMAX**2/DBLE(SCALE))
         XMINV = DBLE(XPR)*(1 + VMIN**2/DBLE(SCALE))
c                  write(6,*) 'pomstr:',XMAX,XMIN
         XR = XMINV*(XMAXV/XMINV)**DBLE(XINT)
         FXT2 = SNGL(XR * DLOG(XMAXV/XMINV))
         IF(XMAXV.GT.XMAX) RETURN
      ENDIF
      T2MN = -(DBLE(SCALE) *(1.D0-XR/DBLE(XPR))
     +         +4.D0*DBLE(ULMASS(211))**2)
      IF(T2MN.LT.T2MAX) THEN
         T2MAX1 = T2MN
      ELSE
         T2MAX1 = T2MAX
      ENDIF

      T2MIN=MP*MP*XR*XR/(1.d0-XR)
      T2 = T2MIN*((T2MAX1/T2MIN)**DBLE(T2INT))
c      write(6,*) ' FXT2 : T2,T2MAX,T2MIN,T2INT ',T2,T2MAX,T2MIN,T2INT
c      T2 = DBLE(T2INT)
      XP2T = DBLE(XPR)
      XX = XP2T/XR

      IF(T2MIN.GE.T2MAX1) THEN
c          write(6,*) ' T2MIN = ',T2MIN,' T2MAX ',T2MAX
         FXT2 = 0.0
         RETURN
      ENDIF
      IF(XP2T.GE.XMAX) THEN
         FXT2 = 0.0
         RETURN
      ENDIF

      if(XR.Le.0.0) THEN
         write(6,*) ' FXT2: xx,XPR,XR ',XX,XPR,XR
         write(6,*) ' FXT2: T2min = ',T2MIN,'T2MAX ',T2MAX
c        FUNX = 0
c        RETURN
      endif
      T2GKI = SNGL(-T2)
      XFGKI = SNGL(XR)
      CALL RAT2DI(KF,XR,-T2,WTDIST)
      CALL RASTFU(KF,SNGL(XX),SCALE,XPQ)

c      write(6,*) ' FXT2: WTDIST =',wtdist
c      write(6,*) 'FXT2 ',xpq(kpart),xpq(-kpart),kpart,xx,scale
      FXT2= FXT2 * SNGL(WTDIST)
      FXT2 = FXT2 * SNGL(T2 * DLOG(T2MAX1/T2MIN))
      IF(IVM.EQ.0) THEN
         WTPI = WTDIST*T2*DLOG(T2MAX1/T2MIN)*XR*DLOG(XMAX/XP2T)
      ELSE
         WTPI = WTDIST*T2*DLOG(T2MAX1/T2MIN)*XR*DLOG(XMAXV/XMINV)
      ENDIF
ctest      WTPI = WTDIST*T2*XR
      IF((NG.GE.30.AND.NG.LT.100).OR.(NG.LT.0.AND.NPOM.LT.0)) THEN
         WGMAX=0.
         DO 10  I=-6,6
            IF(I.EQ.0) GOTO 10
            WGMAX = MAX(WGMAX,XPQ(I))
   10    CONTINUE
         WTPI = WTPI*DBLE(WGMAX)
         IF(WTPI.GT.POM) POM=WTPI
      ELSE
         IF(WTPI.GT.PIM) PIM=WTPI
      ENDIF
      FXT2 = FXT2 * XPQ(KPART)

c      write(6,*) ' FXT2: FXT2 = ',FXT2
      RETURN
      END
