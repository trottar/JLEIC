      SUBROUTINE RAT2IN(KF,XR,XP2,T2MIN,T2,WTDIST)
	Implicit None
c give t distribution in proton
C KF = 100 for pomeron
C KF = 211 for pi
C XR = fractional energy of pomeron
C XP2 = fractional energy of parton to hard interaction
C T2 = t = (p - p')**2 < 0
C WTDIST = weight of distribution
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.

*KEEP,RGPARA.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE


*KEND.
      INTEGER KF
      DOUBLE PRECISION T2,XR,XP2,T2MIN,WTDIST,MP,XE,XMAX,T2CH,beta02
      LOGICAL FIRST

      DATA MP/0.938D0/,FIRST/.TRUE./
      DATA beta02/58.74D0/
      XMAX=1.d0-XF

      WTDIST = 0.D0
c      write(6,*) 'RAT2IN: XR,XP2,T2MIN,T2',XR,XP2,T2MIN,T2
      IF(NPOM.EQ.0) THEN
         IF(FIRST) THEN
            write(6,*) ' streng pomeron is used'
            FIRST = .FALSE.
         ENDIF
c pomeron distribution from streng (hera proc. 1987)
         XE=2.D0*ALPHP*DLOG(1.d0/XR)+RN2
         T2CH = T2MAX
         IF(XE*T2MIN.GT.170.D0) THEN
c                      WRITE(6,*) 'XE..',XE,DEXP(-XE)
            GOTO 10
         ENDIF
         IF(XE*T2MAX.GT.170.D0) T2CH = 100.D0/XE
c         write(6,*) ' XE ',xe,' XR ',XR
c         write(6,*) 'T2MAX ',T2MAX,' T2MIN ',T2MIN
         WTDIST = DEXP(-XE*T2MIN)-DEXP(-XE*T2CH)
c         write(6,*) WTDIST
         WTDIST = WTDIST/XE
         WTDIST = WTDIST*beta02/XR**(1.D0 + 2.D0 * EPSP) /16.d0/PI
         IF(DABS(WTDIST).LT.1.D-20) GOTO 10
      ELSEIF(NPOM.EQ.1) THEN
         IF(FIRST) THEN
            write(6,*) ' Ingelman pomeron is used'
            FIRST = .FALSE.
         ENDIF
c test with Ingelman f_pom distribution
         T2CH = T2MAX
         IF(8.D0*T2MAX.GT.170.D0) T2CH = 170.D0/8.D0
         WTDIST = 6.38D0/8.D0*(DEXP(-8.D0*T2MIN) - DEXP(-8.D0*T2CH))
         WTDIST = WTDIST+ 0.424D0/8.D0*(DEXP(-3*T2MIN) - DEXP(-3*T2CH))
         WTDIST = WTDIST/XR/2.3D0
c end Ingelman
      ELSE
         write(6,*) ' RAT2IN: pomeron distribution ',NPOM,
     +       ' not implemented'
         write(6,*) ' RAT2IN: program stops '
         stop
      ENDIF
      RETURN
   10 WTDIST = 0.D0
      RETURN
      END
