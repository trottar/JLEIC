*CMZ :  2.08/04 22/12/99  15.39.23  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  20.32.53  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE FRAG
	IMPLICIT NONE
*KEEP,RGLUCO.
      REAL PLEPIN,PPIN
      INTEGER KE,KP,KEB,KPH,KGL,KPA,NFRAG,ILEPTO,IFPS,IHF,IALMKT
      INTEGER INTER,ISEMIH
      INTEGER NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT,NFLAV,NFLQCDC
      COMMON/LUCO  /KE,KP,KEB,KPH,KGL,KPA,NFLAV,NFLQCDC
      COMMON/INPU  /PLEPIN,PPIN,NFRAG,ILEPTO,IFPS,IHF,IALMKT,INTER,
     +              ISEMIH
      COMMON/HARD/ NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT
      INTEGER IHFLA
      COMMON/HFLAV/ IHFLA
C      SAVE

*KEEP,RGLUJETS.
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      REAL ULMASS
      DOUBLE PRECISION DOT1,DPLU,DLANGL
      EXTERNAL DLANGL,ULMASS,DPLU,LUCHGE,DOT1
C      SAVE

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
	Integer IJOIN1,IJOIN2,NJOIN,NFLSP,NFLCH,KFL,KFLSP,KINTER
	Integer KFLCH
	REAL PHISPL,PTSPL,CHI,draprn
      DIMENSION IJOIN1(10),IJOIN2(10)
      DATA PHISPL/0./,PTSPL/0./
      DATA NFLSP/0/,NFLCH/0/
      IF(NFRAG.EQ.10.AND.NPOM.NE.20.AND.NPOM.NE.21
     &                .AND.IDIR.EQ.0) CALL PRODIFF
      IF(NPOM.EQ.20.OR.NPOM.EQ.21) THEN
         KFL = KINT(2,2)
         KINTER = K(NIA2,2)
         CALL PYSPLI(KFL,KINTER,KFLCH,KFLSP)
c         write(6,*) ' KFLCH =',KFLCH,' KFLSP =',KFLSP
         K(NF1+1,1)=21
C-- new test
         CHI = draprn()
         IF(IABS(KFLCH).GT.0.AND.IABS(KFLSP).GT.0) THEN
            K(N+1,1)=2
            K(N+1,2)=KFLCH
            P(N+1,1)=P(NF1+1,1)*DBLE((1.0-CHI)+PTSPL*COS(PHISPL))
            P(N+1,2)=P(NF1+1,2)*DBLE((1.0-CHI)+PTSPL*SIN(PHISPL))
            P(N+1,5)=DBLE(ULMASS(KFLCH))
            NFLCH = N+1
            K(N+2,1)=2
            K(N+2,2)=KFLSP
            P(N+2,5)=DBLE(ULMASS(KFLSP))
            P(N+2,1)=P(NF1+1,1)*DBLE(CHI)-DBLE(PTSPL*COS(PHISPL))
            P(N+2,2)=P(NF1+1,2)*DBLE(CHI)-DBLE(PTSPL*SIN(PHISPL))
            NFLSP = N+2
            P(N+1,4)=DBLE((1.0-CHI))*P(NF1+1,4)
            P(N+1,3)=DBLE((1.0-CHI))*P(NF1+1,3)
            P(N+2,4)=DBLE(CHI)*P(NF1+1,4)
            P(N+2,3)=DBLE(CHI)*P(NF1+1,3)
            N = N+2
            NJOIN = 2
            IF(IABS(KFLCH).LT.10.AND.IABS(KFLSP).LT.10) THEN
               IF(KFLCH.GT.0.AND.IABS(KFLSP).GT.0) THEN
                  IJOIN1(1) = NF2
                  IJOIN1(2) = NFLCH
                  IJOIN2(1) = NF1
                  IJOIN2(2) = NFLSP
               ELSE
                  IJOIN1(1) = NF1
                  IJOIN1(2) = NFLCH
                  IJOIN2(1) = NF2
                  IJOIN2(2) = NFLSP
               ENDIF
               CALL LUJOIN(NJOIN,IJOIN1)
               CALL LUJOIN(NJOIN,IJOIN2)
            ELSEIF(IABS(KFLCH).GT.10.AND.IABS(KFLSP).LT.10) THEN
               IJOIN1(1) = NF1
               IJOIN1(2) = NFLSP
               NJOIN = 2
               CALL LUJOIN(NJOIN,IJOIN1)
            ENDIF
         ELSEIF(IABS(KFLCH).EQ.0.AND.IABS(KFLSP).LT.10) THEN
            K(N+1,1)=2
            K(N+1,2)=KFLSP
            P(N+1,1)=P(NF1+1,1)
            P(N+1,2)=P(NF1+1,2)
            P(N+1,3)=P(NF1+1,3)
            P(N+1,4)=P(NF1+1,4)
            P(N+1,5)=DBLE(ULMASS(KFLSP))
            NFLSP = N+1
            NJOIN = 2
            IJOIN1(1) = NF1
            IJOIN1(2) = NFLSP
            N=N+1
            CALL LUJOIN(NJOIN,IJOIN1)
         ELSEIF(IABS(KFLCH).LT.10.AND.IABS(KFLSP).EQ.0) THEN
            K(N+1,1)=2
            K(N+1,2)=KFLCH
            P(N+1,1)=P(NF1+1,1)
            P(N+1,2)=P(NF1+1,2)
            P(N+1,3)=P(NF1+1,3)
            P(N+1,4)=P(NF1+1,4)
            P(N+1,5)=DBLE(ULMASS(KFLCH))
            NFLCH = N+1
            NJOIN = 2
            IJOIN1(1) = NF1
            IJOIN1(2) = NFLCH
            N=N+1
            CALL LUJOIN(NJOIN,IJOIN1)
         ELSEIF(IABS(KFLCH).GT.10.AND.IABS(KFLSP).LT.10) THEN
            NJOIN = 3
            IJOIN1(1) = NF1
            IJOIN1(2) = NFLSP
            IJOIN1(3) = NF2
            CALL LUJOIN(NJOIN,IJOIN1)
         ELSEIF(IABS(KFLCH).LT.10.AND.IABS(KFLSP).GT.10) THEN
            NJOIN = 3
            IJOIN1(1) = NF1
            IJOIN1(2) = NFLCH
            IJOIN1(3) = NF2
            CALL LUJOIN(NJOIN,IJOIN1)
         ELSE
            write(6,*) ' this configuration not implemented'
         ENDIF
      ENDIF
      CALL DUPREP(NF1)
      RETURN
      END
