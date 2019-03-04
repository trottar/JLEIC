*CMZ :  2.08/05 30/03/2000  16.49.31  by  Hannes Jung
*CMZ :  2.08/02 04/11/99  16.39.09  by  Hannes Jung
*-- Author :    Hannes Jung   08/09/99
      subroutine rghaprep
      implicit none
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

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
*KEND.
      Real     remPARJ32
      Real     remPARJ(5)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
c henry's commons
      Double Precision  Stot, Q2_henry, W2, Mx2, AJAC
      COMMON /CQ2W2MX/  Stot, Q2_henry, W2, Mx2, AJAC
      real rgmx2
      IF(FIRST) THEN
         first = .FALSE.
         remPARJ32 = PARJ(32)
      ELSE
         PARJ(32)=remPARJ32
      ENDIF
      if(ipro.eq.30) then
         rgmx2 = real(mx2)
      else
         rgmx2 = shh
      endif
      If (RGMX2.LT.4.0.AND.IDIR.EQ.0) Then
         If(ABS(KPA).LE.2) Then
            PARJ(32) = 0.35
            GOTO 190
         ElseIf(ABS(KPA).EQ.3) Then
            PARJ(32) = 0.1
            GOTO 190
         EndIf
  190    Continue
*
         remPARJ(1) = PARJ(11)
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         remPARJ(4) = PARJ(27)
         remPARJ(5) = PARJ(28)
*
         PARJ(11) = 1.0
         PARJ(12) = 1.0
         PARJ(13) = 1.0
         PARJ(27) = 0.1
         PARJ(28) = 0.9
*
         CALL DUPREP(0)
         PARJ(11) = remPARJ(1)
         PARJ(12) = remPARJ(2)
         PARJ(13) = remPARJ(3)
         PARJ(27) = remPARJ(4)
         PARJ(28) = remPARJ(5)
c      ElseIf (RGMX2.GT.9.0.AND.RGMX2.LT.14.AND. ABS(KPA) .EQ.4.AND.IDIR
c     +.EQ.0) Then
      ElseIf (RGMX2.GT.9.0.AND. ABS(KPA) .EQ.4.AND.IDIR
     +.EQ.0) Then
         remPARJ(2) = PARJ(12)
         remPARJ(3) = PARJ(13)
         PARJ(13) = 1.0
         CALL DUPREP(0)
         PARJ(13) = remPARJ(3)
         PARJ(12) = remPARJ(2)

      Else

         CALL DUPREP(0)


      ENDIF
      return
      end
