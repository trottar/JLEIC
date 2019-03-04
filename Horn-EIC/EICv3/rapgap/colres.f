*CMZ :  2.07/03 15/05/99  14.31.08  by  Hannes Jung
*-- Author :    Hannes Jung   05/03/97
      SUBROUTINE COLRES

      IMPLICIT NONE
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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
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

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


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
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      INTEGER LST,IRES
      COMMON/EPPARA/LST(30),IRES(2)
c..hju
      INTEGER ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG
      COMMON /COLCON/ ICOLORA,IRESPRO,IRPA,IRPB,IRPC,IRPD,IRPE,IRPF,IRPG

      INTEGER KS,NS
      COMMON /COLR/ KS(40,5),NS


      K(NS+7,2)=21
ccc      write(6,*) ' colres IRESPRO ',irespro
      IF(IRESPRO.EQ.1) THEN
         K(NS+3,4)=21
         K(NS+3,5)=27
c qqbar event
         IF(K(NS+5,2).GT.0) THEN
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=0
            K(NS+6,4)=0
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSE
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=0
            K(NS+5,4)=0
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ENDIF
      ELSEIF(IRESPRO.EQ.2) THEN
c glu glu event

         IF(ICOLORA.EQ.1) THEN
c color flow A
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+7)*MSTU(5)
            K(NS+5,5)=(NS+6)*MSTU(5)
            K(NS+6,4)=(NS+5)*MSTU(5)
            K(NS+6,5)=(NS+7)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+25
            K(NS+7,5)=(NS+3)*MSTU(5)+26
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            K(NS+1,4) = 27
            K(NS+1,5) = 23
            K(NS+3,4)=21
            K(NS+3,5)=27
            K(NS+5,4)=(NS+6)*MSTU(5)
            K(NS+5,5)=(NS+7)*MSTU(5)
            K(NS+6,4)=(NS+7)*MSTU(5)
            K(NS+6,5)=(NS+5)*MSTU(5)
            K(NS+7,4)=(NS+1)*MSTU(5)+26
            K(NS+7,5)=(NS+3)*MSTU(5)+25
         ELSEIF(ICOLORA.EQ.3) THEN
c color flow C
            K(NS+1,4) = 26
            K(NS+1,5) = 25
            K(NS+3,4)=25
            K(NS+3,5)=26
            K(NS+5,4)=(NS+3)*MSTU(5)
            K(NS+5,5)=(NS+1)*MSTU(5)
            K(NS+6,4)=(NS+1)*MSTU(5)
            K(NS+6,5)=(NS+3)*MSTU(5)
            K(NS+7,4)=0
            K(NS+7,5)=0
            K(NS+7,1)=0
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF

      ELSEIF(IRESPRO.EQ.3) THEN
c qg event
C color configuartion A
         IF(ICOLORA.EQ.1) THEN
C antq from photon
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 23
                  K(NS+3,4)=21
                  K(NS+3,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from photon
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 23
                  K(NS+1,5) = 0
                  K(NS+3,4)=27
                  K(NS+3,5)=21
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 21
                  K(NS+1,4)=23
                  K(NS+1,5)=27
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+5,4)=(NS+6)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=(NS+5)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 21
                  K(NS+3,5) = 0
                  K(NS+1,4)=27
                  K(NS+1,5)=23
                  K(NS+6,4)=(NS+5)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=(NS+6)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF
C color configuartion B
         ELSE
C antq from photon
            IF(K(21,2).LT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 26
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+1,4) = 0
                  K(NS+1,5) = 25
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,2) = K(NS+6,2)
                  K(NS+7,1)=0
c               write(6,*) ' 3rd '
               ENDIF
c  q from photon
            ELSEIF(K(21,2).GT.0.AND.K(23,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+1,4) = 26
                  K(NS+1,5) = 0
                  K(NS+3,4)=25
                  K(NS+3,5)=26
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+1,4) = 25
                  K(NS+1,5) = 0
                  K(NS+3,4)=26
                  K(NS+3,5)=25
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
C antq from proton
            ELSEIF(K(23,2).LT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).LT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 26
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 1st '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).LT.0) THEN
                  K(NS+3,4) = 0
                  K(NS+3,5) = 25
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=(NS+3)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 3rd '
               ENDIF
c  q from proton
            ELSEIF(K(23,2).GT.0.AND.K(21,2).EQ.21) THEN
               IF(K(25,2).GT.0.AND.K(26,2).EQ.21) THEN
                  K(NS+3,4) = 26
                  K(NS+3,5) = 0
                  K(NS+1,4)=25
                  K(NS+1,5)=26
                  K(NS+5,4)=(NS+1)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=(NS+3)*MSTU(5)
                  K(NS+6,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+5,2)
c               write(6,*) ' 2nd '
               ELSEIF(K(25,2).EQ.21.AND.K(26,2).GT.0) THEN
                  K(NS+3,4) = 25
                  K(NS+3,5) = 0
                  K(NS+1,4)=26
                  K(NS+1,5)=25
                  K(NS+6,4)=(NS+1)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+5,4)=(NS+3)*MSTU(5)
                  K(NS+5,5)=(NS+1)*MSTU(5)
                  K(NS+7,4)=0
                  K(NS+7,5)=0
                  K(NS+7,1)=0
                  K(NS+7,2) = K(NS+6,2)
c               write(6,*) ' 4th '
               ENDIF
            ENDIF

         ENDIF
ccc         CALL DULIST(2)
      ELSEIF(IRESPRO.EQ.4) THEN
c qqbar --> g g
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+26
               K(NS+7,5)=(NS+1)*MSTU(5)+25
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+25
               K(NS+7,5)=(NS+3)*MSTU(5)+26
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               K(NS+5,4)=(NS+7)*MSTU(5)
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=(NS+7)*MSTU(5)
               K(NS+7,4)=(NS+3)*MSTU(5)+25
               K(NS+7,5)=(NS+1)*MSTU(5)+26
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=(NS+7)*MSTU(5)
               K(NS+6,4)=(NS+7)*MSTU(5)
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=(NS+1)*MSTU(5)+26
               K(NS+7,5)=(NS+3)*MSTU(5)+25
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.5) THEN
c qqbar --> q qbar
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 23
               K(NS+3,4)=21
               K(NS+3,5)=0
               K(NS+5,4)=0
               K(NS+5,5)=(NS+6)*MSTU(5)
               K(NS+6,4)=(NS+5)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 23
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=21
               K(NS+5,4)=(NS+6)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=0
               K(NS+6,5)=(NS+5)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 27
               K(NS+3,4)=27
               K(NS+3,5)=0
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+3)*MSTU(5)+25
                  K(NS+7,5)=(NS+1)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+3)*MSTU(5)+26
                  K(NS+7,5)=(NS+1)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 27
               K(NS+1,5) = 0
               K(NS+3,4)=0
               K(NS+3,5)=27
               IF(K(25,2).GT.0) THEN
                  K(NS+5,4)=(NS+7)*MSTU(5)
                  K(NS+5,5)=0
                  K(NS+6,4)=0
                  K(NS+6,5)=(NS+7)*MSTU(5)
                  K(NS+7,4)=(NS+1)*MSTU(5)+25
                  K(NS+7,5)=(NS+3)*MSTU(5)+26
               ELSE
                  K(NS+5,4)=0
                  K(NS+5,5)=(NS+7)*MSTU(5)
                  K(NS+6,4)=(NS+7)*MSTU(5)
                  K(NS+6,5)=0
                  K(NS+7,4)=(NS+1)*MSTU(5)+26
                  K(NS+7,5)=(NS+3)*MSTU(5)+25
               ENDIF
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.6) THEN
c q q --> q q
         IF(ICOLORA.EQ.1) THEN
c color flow A
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 26
               K(NS+1,5) = 0
               K(NS+3,4)=25
               K(NS+3,5)=0
               K(NS+5,4)=(NS+3)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+1)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 26

               K(NS+3,4)=0
               K(NS+3,5)=25
               K(NS+5,4)=0
               K(NS+5,5)=(NS+3)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+1)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow A, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSEIF(ICOLORA.EQ.2) THEN
c color flow B
            IF(K(23,2).GT.0) THEN
               K(NS+1,4) = 25
               K(NS+1,5) = 0

               K(NS+3,4)=26
               K(NS+3,5)=0
               K(NS+5,4)=(NS+1)*MSTU(5)
               K(NS+5,5)=0
               K(NS+6,4)=(NS+3)*MSTU(5)
               K(NS+6,5)=0
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, quark'
            ELSEIF(K(23,2).LT.0) THEN
               K(NS+1,4) = 0
               K(NS+1,5) = 25

               K(NS+3,4)=0
               K(NS+3,5)=26
               K(NS+5,4)=0
               K(NS+5,5)=(NS+1)*MSTU(5)
               K(NS+6,4)=0
               K(NS+6,5)=(NS+3)*MSTU(5)
               K(NS+7,4)=0
               K(NS+7,5)=0
               K(NS+7,1)=0
c            write(6,*) ' flow B, antiquark'
            ELSE
               write(6,*) ' wrong flavor',K(23,2),' for IRESPRO=',
     +         IRESPRO
               write(6,*) ' color conf: ',ICOLORA
            ENDIF
         ELSE
            write(6,*) ' wrong color configuartion ',ICOLORA
         ENDIF
      ELSEIF(IRESPRO.EQ.7) THEN
c q q --> q q color singlet exchange
         K(NS+1,4) = 25
         K(NS+1,5) = 25
         K(NS+3,4)=26
         K(NS+3,5)=26
         K(NS+5,4)=(NS+1)*MSTU(5)
         K(NS+5,5)=(NS+1)*MSTU(5)
         K(NS+6,4)=(NS+3)*MSTU(5)
         K(NS+6,5)=(NS+3)*MSTU(5)
         K(NS+7,4)=0
         K(NS+7,5)=0
         K(NS+7,1)=0
      ELSE
         write(6,*) ' colres : irespro not implemented ',irespro
      ENDIF
      RETURN
      END
