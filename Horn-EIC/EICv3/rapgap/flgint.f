*CMZ :  2.06/51 10/10/98  10.30.09  by  Hannes Jung
*CMZ :  2.06/44 28/07/98  18.49.16  by  Hannes Jung
*CMZ :  2.06/39 16/07/98  11.40.09  by  Hannes Jung
*CMZ :  2.06/38 05/07/98  21.25.53  by  Hannes Jung
*CMZ :  6.22/00 20/07/95  17.30.36  by  Martin Hampel
*-- Author :
C **********************************************************************

      FUNCTION FLGINT(Z)

C...Gluon contribution integrand to QCD longitudinal structure function.
      DOUBLE PRECISION PARL
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
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

*KEND.
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      IF(IHF.EQ.0) THEN
         CALL LNSTRF(Z,Q2,XPQ)
         FLGINT=20./9.*PARL(25)/PI*(X/Z)**2*(1.-X/Z)/Z*XPQ(0)
         NPASS=NPASS+1
      ELSEIF(IHF.EQ.1) THEN
         FLGINT=FLHF(Z,X,Q2)
         NPASS=NPASS+1
      ELSE
         write(6,*) ' FLGINT wrong IHF = ',IHF
         STOP
      ENDIF


      RETURN
      END
