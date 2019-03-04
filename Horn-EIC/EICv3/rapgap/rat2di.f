*CMZ :  2.08/05 10/03/2000  12.02.57  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  09.07.52  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  10.08.10  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE RAT2DI(KF,XR,T2,WTDIST)
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
      DOUBLE PRECISION T2,XR,WTDIST,DELTA2,MP,XE,BETA02,WMOPOL,WDIPOL
      DOUBLE PRECISION WEXPOL,WEXTHO,WEXPHO
      DOUBLE PRECISION GPPPI,mpi,ALAMD,RNHOL
      LOGICAL FIRST
      DATA DELTA2/3.26D0/,MP/0.938D0/,FIRST/.TRUE./
      DATA GPPPI/13.5/,mpi/0.139/,ALAMD/0.039D0/,RNHOL/0.93/
      DATA beta02/58.74D0/
      WTDIST = 0.D0
c        write(6,*) ' RAT2DI : NG,NPOM,KF,IPRO',NG,NPOM,KF,IPRO
c      IF(IWEI.EQ.1) THEN
c        write(6,*) ' RAT2DI : NG,NPOM,KF,IPRO',NG,NPOM,KF,IPRO
c        ENDIF
      IF(NPOM.EQ.0) THEN
         IF(FIRST) THEN
            write(6,*) ' streng pomeron is used'
            FIRST = .FALSE.
         ENDIF
c pomeron distribution from streng (hera proc. 1987)
         XE=RN2*DABS(T2)
         IF(XE.GT.170.D0) THEN
C                      WRITE(6,*) 'XE..',XE,DEXP(-XE)
            GOTO 10
         ENDIF
         WTDIST = DEXP(-XE)
         IF(WTDIST.LT.1.D-50) GOTO 10
         WTDIST = WTDIST*(XR**(1.D0-2.D0*(1.D0 + EPSP + ALPHP*T2)))
         WTDIST = WTDIST*beta02/16.D0/PI
c         write(6,*) 'rat2di : NPOM=',NPOM
c          write(6,*) 'rat2di: ',epsp,alphp,t2,rn2
      ELSEIF(NPOM.EQ.1) THEN
         IF(FIRST) THEN
            write(6,*) ' Ingelman pomeron is used'
            FIRST = .FALSE.
         ENDIF
c test with Ingelman f_pom distribution
         WTDIST =(6.38D0*DEXP(8.D0*T2) + 0.424D0*DEXP(3.D0*T2))/XR/2.3D0
c         write(6,*) 'rat2di : NPOM=',NPOM
c end Ingelman
      ELSEIF(NPOM.EQ.2) THEN
         IF(FIRST) THEN
            write(6,*) ' Donnachie Landshoff pomeron is used'
            FIRST = .FALSE.
         ENDIF
c test with Donnachie Landshoff f_pom distribution
         WTDIST = 9.D0*DELTA2/4.D0/PI**2
         WTDIST = WTDIST*(4.D0*MP**2 - 2.8D0*T2)/(4.D0*MP**2 - T2)
         WTDIST = WTDIST/(1.D0 - T2/0.7D0)
         WTDIST = WTDIST*(XR**(1.D0-2.D0*(1.D0 + EPSP + ALPHP*T2)))
c         write(6,*) 'rat2di : NPOM=',NPOM
      ELSEIF(NPOM.EQ.20.OR.NPOM.EQ.21) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' pi exchange is used '
            FIRST = .FALSE.
         ENDIF
changed according to B.List/Holtmann
         WTDIST = -2.D0*GPPPI**2 *XR/16.D0/PI**2 *T2/(T2 - mpi**2)**2
c different choices for pion t dependence a la Frankfurt Mankiewiecz Strikman
c  Z.Phys. A 334 (1989) 343
c dipole form factor
         WDIPOL = (1.D0 - mpi**2/0.9D0/0.9D0)**2
     +   /(1.D0 -T2/0.9D0/0.9D0)**2
         WDIPOL = WDIPOL **2
c monopol form factor
         WMOPOL = (1.D0 - mpi**2/0.5D0/0.5D0)/(1.D0 -T2/0.5D0/0.5D0)
         WMOPOL = WMOPOL **2
c exp. form factor
         WEXPOL = dexp(1.8D0 * mpi**2) * dexp(1.8D0 *T2)
         WEXPOL = WEXPOL **2
c exp form factor from Thomas Phys. Lett. B 126 (1983) 97
         WEXTHO =  DABS(DEXP(-2.D0*ALAMD *(-T2 + mpi**2)/mpi**2))
c exp form factor  from Holtmann et al. Phys. Lett. B338 (1994) 363
         WEXPHO = DEXP(RNHOL**2*(T2-mpi**2)/XR)
c	   write(6,*) ' rat2di ',wtdist,wexpho,xr,t2,rnhol,mpi
c.....................................................................
         WTDIST = WTDIST * WEXPHO
c         write(6,*) 'rat2di : NPOM=',NPOM
         IF(NPOM.EQ.21) THEN
c construct t distribution for pi0 exchange
            WTDIST = WTDIST/2.D0
         ENDIF
      ELSEIF(NPOM.EQ.30) THEN
         IF(FIRST) THEN
            write(6,*) ' RAT2DI: Nikolaev Zakharov process selected'
            FIRST = .FALSE.
         ENDIF
         WTDIST = 1.0D0
      ELSEIF(NPOM.EQ.40) THEN
         IF(FIRST) THEN
            write(6,*) ' RAT2DI: M. Wuesthoff model used'
            FIRST = .FALSE.
         ENDIF
         WTDIST = 1.0D0
      ELSEIF(NPOM.EQ.41) THEN
         IF(FIRST) THEN
            write(6,*) ' RAT2DI: Bartels Lotter Wuesthoff calc. used'
            FIRST = .FALSE.
         ENDIF
         WTDIST = 1.0D0
      ELSEIF(NPOM.EQ.42) THEN
         IF(FIRST) THEN
            write(6,*) ' RAT2DI: M. Diehl calc. used'
            FIRST = .FALSE.
         ENDIF
         WTDIST = 1.0D0
      ELSEIF(NPOM.EQ.45) THEN
         IF(FIRST) THEN
            WRITE(6,*) ' RASTFU: Buchmueller/McDermott/Hebecker
     &        calc. used '
            FIRST = .FALSE.
         ENDIF
         WTDIST = 1.0D0
      ELSEIF(NPOM.LT.0) THEN
c user supplied pomeron distribution
         WTDIST = 1.D0
c         write(6,*) 'rat2di : NPOM=',NPOM
      ELSE
         write(6,*) ' RAT2DI: pomeron distribution ',NPOM,
     +              ' not implemented'
         write(6,*) ' RAT2DI: program stops '
         stop
      ENDIF
c	write(6,*) ' RAT2DI: ',xr,t2,wtdist
      RETURN
   10 WTDIST = 0.D0
      RETURN
      END
