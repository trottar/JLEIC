*CMZ :  2.08/00 06/06/99  15.45.42  by  Hannes Jung
*CMZ :  2.06/51 10/10/98  10.29.55  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.48  by  Hannes Jung
*CMZ :  2.06/21 28/11/97  11.56.06  by  Hannes Jung
*CMZ :  2.06/08 10/10/97  17.15.19  by  Hannes Jung
*CMZ :  2.06/00 17/07/97  12.49.00  by  Hannes Jung
*CMZ :  2.04/00 10/12/96  11.15.18  by  Hannes Jung
*CMZ :  2.01/14 03/05/96  10.11.16  by  Hannes Jung
*CMZ :  2.01/13 28/04/96  16.57.17  by  Hannes Jung
*CMZ :  2.01/09 10/02/96  19.20.34  by  Hannes Jung
*CMZ :  2.01/00 01/01/96  17.20.46  by  Hannes Jung
*CMZ :  2.00/30 29/12/95  14.46.41  by  Hannes Jung
*CMZ :  2.00/28 25/11/95  17.12.18  by  Hannes Jung
*CMZ :  2.00/27 22/11/95  16.24.57  by  Hannes Jung
*CMZ :  2.00/25 20/11/95  16.31.16  by  Hannes Jung
*CMZ :  2.00/23 19/11/95  21.32.59  by  Hannes Jung
*CMZ :  1.04/00 11/01/95  14.35.16  by  Hannes Jung
*CMZ :  1.03/01 03/04/94  16.45.58  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE AREXEC
C...ARiadne subroutine EXECute ariadne

C...The Main driver routine in Ariadne.


      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT DOUBLE PRECISION (B)
      IMPLICIT LOGICAL (Q)
      PARAMETER(MAXDIP=500,MAXPAR=500,MAXSTR=100)
      COMMON /ARPART/ BP(MAXPAR,5),IFL(MAXPAR),IEX(MAXPAR),QQ(MAXPAR),
     +                IDI(MAXPAR),IDO(MAXPAR),INO(MAXPAR),IPART
      SAVE /ARPART/

      COMMON /ARDIPS/ BX1(MAXDIP),BX3(MAXDIP),PT2IN(MAXDIP),
     +                SDIP(MAXDIP),IP1(MAXDIP),IP3(MAXDIP),
     +                AEX1(MAXDIP),AEX3(MAXDIP),QDONE(MAXDIP),
     +                QEM(MAXDIP),IRAD(MAXDIP),ISTR(MAXDIP),IDIPS
      SAVE /ARDIPS/

      COMMON /ARSTRS/ IPF(MAXSTR),IPL(MAXSTR),IFLOW(MAXSTR),
     +                PT2LST,IMF,IML,IO,QDUMP,ISTRS
      SAVE /ARSTRS/

      COMMON /ARDAT1/ PARA(40),MSTA(40)
      SAVE /ARDAT1/

      COMMON /ARDAT2/ PQMAS(10)
      SAVE /ARDAT2/

      COMMON /ARDAT3/ IWRN(40)
      SAVE /ARDAT3/

      COMMON /ARINT2/ DBEX,DBEY,DBEZ,PHI,THE
      SAVE /ARINT2/

*KEEP,RGLUPARM.
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


*KEND.
      REAL P,V
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      DOUBLE PRECISION PARL
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,XQ2,U
      SAVE /RAPTOU/
*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEND.

      COMMON /PYSUBS/ MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
      SAVE /PYSUBS/

      COMMON /PYINT1/ MINT(400),VINT(400)
      SAVE /PYINT1/
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

*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEEP,RGPARAS.
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

*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEND.
C...Step counter
      MSTA(4)=MSTA(4)+1

C...Reset error log
      MSTA(13)=0

C...Error if ARINIT has not been called
      IF(MSTA(2).EQ.0) CALL ARERRM('AREXEC',12,0)

C...If ariadne mode just pass event through to ARPARS
      IF(MSTA(1).EQ.0) THEN
         CALL ARPARS(1,N)

C...If JETSET mode should work by just passing event on to ARPARS
      ELSEIF(MSTA(1).EQ.1) THEN
         CALL ARPARS(1,N)

C...If PYTHIA mode tag extended partons etc.
      ELSEIF(MSTA(1).EQ.2) THEN

         ISUB=MINT(1)
         IF(ISUB.NE.11.AND.ISUB.NE.12.AND.ISUB.NE.13.AND.ISUB.NE.28.AND
     +   .ISUB.NE.53.AND.ISUB.NE.68) CALL ARERRM('AREXEC',20,0)

         IFIRST=1
         ILAST=N

         DO 10 I=IFIRST,ILAST
            IF(K(I,1).GT.2) GOTO 10
            CALL ARGTYP(I,ITYP)
            IF(ITYP.EQ.0) GOTO 10
            IF(K(I,3).EQ.1.OR.K(I,3).EQ.2) THEN
               K(I,4)=1
            ELSE
               K(I,4)=0
            ENDIF
   10    CONTINUE

         CALL ARPARS(IFIRST,ILAST)

C...If LEPTO mode tag extended partons
      ELSEIF(MSTA(1).EQ.3) THEN
         IF(LST(24).EQ.1) THEN

C...Boost to hadronic cm to avoid precision problems
            DEL=DBLE(P(5,4)+P(6,4))
            DBXL=DBLE(P(5,1)+P(6,1))/DEL
            DBYL=DBLE(P(5,2)+P(6,2))/DEL
            DBZL=DBLE(P(5,3)+P(6,3))/DEL
            CALL DUDBRB(5,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)

            IF(MSTA(30).LT.2) THEN
               K(5,4)=0
            ELSE
               K(5,4)=3
               PARA(13)=SQRT(XQ2)
            ENDIF
            IF(MSTA(30).EQ.0) THEN
               K(6,4)=1
            ELSE
               K(6,4)=2
               PARA(12)=PARA(11)/(1.0-X)
            ENDIF
            CALL ARPARS(5,6)
            CALL DUDBRB(5,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ELSEIF(LST(24).EQ.3) THEN

C...Boost to hadronic cm to avoid precision problems
            DEL=DBLE(P(5,4)+P(6,4)+P(7,4)+P(8,4))
            DBXL=DBLE(P(5,1)+P(6,1)+ P(7,1)+P(8,1))
     +      /DEL
            DBYL=DBLE(P(5,2)+P(6,2)+ P(7,2)+P(8,2))
     +      /DEL
            DBZL=(DBLE(P(5,3))+DBLE(P(6,3))+ DBLE(P(7,3))+DBLE(P(8,4)))
     +      /DEL
            CALL DUDBRB(5,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)

            IF(MSTA(30).LT.2) THEN
               K(5,4)=0
            ELSE
               K(5,4)=3
               PARA(13)=SQRT(XQ2)
            ENDIF
            IF(MSTA(30).EQ.0) THEN
               K(6,4)=1
            ELSE
               K(6,4)=2
               PARA(12)=PARA(11)/(1.0-X)
            ENDIF
            CALL ARPARS(5,6)
            IF(MSTA(30).LT.2) THEN
               K(7,4)=0
            ELSE
               K(7,4)=3
               PARA(13)=SQRT(XQ2)
            ENDIF
            IF(MSTA(30).EQ.0) THEN
               K(8,4)=1
            ELSE
               K(8,4)=2
               PARA(12)=PARA(11)/(1.0-X)
            ENDIF
            CALL ARPARS(7,8)
            CALL DUDBRB(5,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ENDIF
      ELSEIF(MSTA(1).EQ.4) THEN
c         write(6,*) ' in arexec ',idir,ipro
         INPOM = 2
         XRA = Q2/dble(yy)/sss
         IF(IPRO.EQ.12) THEN
c here QPM
            N1 = 0
            N2 = 0
            DO 30 I=3,N
               IF(IDIR.EQ.0.AND.IABS(K(I,2)).EQ.100) INPOM = I
               IF(K(I,1).EQ.2.AND.N1.EQ.0) THEN
                  N1 = I
                  DO 20 J=N1,N
                     IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N2.EQ.0)
     +               THEN
                        N2 = J+1
                        GOTO 30
                     ENDIF
   20             CONTINUE
               ENDIF
   30       CONTINUE
C...Boost to hadronic cm to avoid precision problems
            DEL=DBLE(P(NIA1,4)+P(INPOM,4))
            DBXL=DBLE(P(NIA1,1)+P(INPOM,1))/DEL
            DBYL=DBLE(P(NIA1,2)+P(INPOM,2))/DEL
            DBZL=DBLE(P(NIA1,3)+P(INPOM,3))/DEL
            CALL DUDBRB(0,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)
            IF(N1.NE.0.AND.N2.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N1,4)=0
               ELSE
                  K(N1,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N2,4)=1
               ELSE
                  K(N2,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
c         write(6,*) ' arexec ipro=12,n1,n2 ',n1,n2
               CALL ARPARS(N1,N2)
            ENDIF
            CALL DUDBRB(0,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ELSEIF(IPRO.EQ.13.OR.IPRO.EQ.14) THEN
c here  BGF
            N1 = 0
            N2 = 0
            DO 50 I=3,N
               IF(IDIR.EQ.0.AND.IABS(K(I,2)).EQ.100) INPOM = I
               IF(K(I,1).EQ.2.AND.N1.EQ.0) THEN
                  N1 = I
                  DO 40 J=N1,N
                     IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N2.EQ.0)
     +               THEN
                        N2 = J + 1
                        GOTO 60
                     ENDIF
   40             CONTINUE
               ENDIF
   50       CONTINUE
   60       CONTINUE
C...Boost to hadronic cm to avoid precision problems
            DEL=DBLE(P(NIA1,4)+P(INPOM,4))
            DBXL=DBLE(P(NIA1,1)+P(INPOM,1))/DEL
            DBYL=DBLE(P(NIA1,2)+P(INPOM,2))/DEL
            DBZL=DBLE(P(NIA1,3)+P(INPOM,3))/DEL
            CALL DUDBRB(0,N,0.D0,0.D0,-DBXL,-DBYL,-DBZL)
            IF(MSTA(30).EQ.0) THEN
               K(NIA2+1,4)=1
            ENDIF
            IF(N2.NE.0) THEN
               NN2 = N2
               N3 = 0
               N4 = 0
               DO 80 I=NN2,N
                  IF(K(I,1).EQ.2.AND.N3.EQ.0) THEN
                     N3 = I
                     DO 70 J=N3,N
                        IF(K(J,1).EQ.2.AND.K(J+1,1).EQ.1.AND.N4.EQ.0)
     +                  THEN
                           N4 = J+1
                           GOTO 90
                        ENDIF
   70                CONTINUE
                  ENDIF
   80          CONTINUE
            ENDIF
   90       CONTINUE
c         write(6,*) ' arexec ipro=13,14  ,n1,n2 ',n1,n2
            IF(N1.NE.0.AND.N2.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N1,4)=0
                  IF(IDIR.EQ.0.AND.(N1+1).LT.N2) K(N2,4)=0
               ELSE
                  K(N1,4)=3
                  IF(IDIR.EQ.0.AND.(N1+1).LT.N2) K(N2,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N1+1,4)=1
               ELSE
                  K(N1+1,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
               CALL ARPARS(N1,N2)
            ENDIF

c         write(6,*) ' arexec 2nd ipro=13,14  ,n3,n4 ',n3,n4
            IF(N3.NE.0.AND.N4.NE.0) THEN
               IF(MSTA(30).LT.2) THEN
                  K(N3,4)=0
               ELSE
                  K(N3,4)=3
                  PARA(13)=SQRT(SNGL(Q2))
               ENDIF
               IF(MSTA(30).EQ.0) THEN
                  K(N4,4)=1
               ELSE
                  K(N4,4)=2
                  PARA(12)=PARA(11)/(1.0-XRA)
                  IF(IDIR.EQ.0) PARA(12)=PARA(11)/(1.0-XRA/XFGKI)
               ENDIF
               CALL ARPARS(N3,N4)
            ENDIF

            CALL DUDBRB(0,N,0.D0,0.D0,DBXL,DBYL,DBZL)
         ENDIF
      ENDIF
C...Perform fragmentation if requested
      IF(MSTA(5).EQ.1) CALL LUEXEC
      RETURN

C**** END OF AREXEC ****************************************************
      END
