*CMZ :  2.08/04 22/12/99  15.39.27  by  Hannes Jung
*CMZ :  2.08/02 14/07/99  18.40.55  by  Hannes Jung
*CMZ :  2.08/01 23/06/99  07.03.33  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  08.44.27  by  Hannes Jung
*-- Author :
      SUBROUTINE QCDMIX(X,IERR)
      IMPLICIT None
      Double Precision X,XV,WMAX
      Integer NDIMEN,LST,IRES,IGENFL,IERR,IMIX
      COMMON /EPPARA/ LST(30),IRES(2)
      COMMON /OALPHAS/ WMAX,IMIX
      COMMON/GENWEI/IGENFL
      DIMENSION X(20)
      COMMON /XVAL/ XV(20),NDIMEN
      Double Precision draprn
	Real RQPM,RQQB,RQQBC,RQQBB,RQCDC,RQQBT,RQPMT
      REAL FL,FU,PQCDQQB,XEPS,PQCDI,PQCD

      Double Precision XG(20)
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
*KEEP,RGQCDGRID.
	Integer NY,NQ,NPPO
      PARAMETER (NY=40,NQ=20,NPPO=20)
cccc      PARAMETER (NY=80,NQ=80)
      DOUBLE PRECISION QY,QQ,QPOM,QPM,QQB,QQBH,QQBB,QCDC
      DOUBLE PRECISION QPMDF,QQBDF,QQBHDF,QQBBDF,QCDCDF
      DOUBLE PRECISION QPMPI,QQBPI,QQBHPI,QQBBPI,QCDCPI
      COMMON/QCDGRI/QY(NY),QQ(NQ),QPOM(NPPO),
     + QPMDF(NY,NQ,NPPO),QQBDF(NY,NQ,NPPO),QQBHDF(NY,NQ,NPPO),
     + QQBBDF(NY,NQ,NPPO),QCDCDF(NY,NQ,NPPO),
     + QPMPI(NY,NQ,NPPO),QQBPI(NY,NQ,NPPO),QQBHPI(NY,NQ,NPPO),
     + QQBBPI(NY,NQ,NPPO),QCDCPI(NY,NQ,NPPO),
     + QPM(NY,NQ),QQB(NY,NQ),QQBH(NY,NQ),QQBB(NY,NQ),QCDC(NY,NQ)
*KEEP,RGRAHER.
      REAL XPQDIF,XPQPI
	Integer IHERPYS
	Integer NBQ2,NBX
      PARAMETER (NBQ2=20)
      PARAMETER (NBX=20)
      COMMON /RAHER/ IHERPYS,XPQDIF(-6:6,NBX,NBQ2),XPQPI(-6:6,NBX,NBQ2)
*KEEP,RGHS45.
      INTEGER IHERAC,KPAHS
      DOUBLE PRECISION YHS,XHS,Q2HS
      COMMON /HS45/ IHERAC,KPAHS,YHS,XHS,Q2HS
*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

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


*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
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

*KEEP,RGFULL.
      INTEGER IFULL,IQCDGRID
      COMMON /OALPINI/ IFULL,IQCDGRID
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.
      Integer NDIM,NPOIN
      COMMON/DIVO/ NDIM,NPOIN
      REAL ULMASS
      REAL SLO,SHI
      REAL SNGL
      Double Precision WMG,WDUM,FXN1,DSMALL,WMAX13,WMAX14,WMAX15
      Double Precision YD,QD,X1P,X2P,WPA,SM,XM_TEST,SMIN,WTEST
      Double Precision PQPM,PQQB,PQQBC,PQQBB,PQCDC,PTEST
      Double Precision X1,Y1,xm1,xm2,xm3,xm4
      DIMENSION WMG(4,NY,NQ)
      Integer K,I,J,NERP,NX1,NX2,NCMAX,IHERPYSO,IDIRO,IGENFLO
      Integer IY,IQ,IPO,IPG,NNX1,NNX2,IFAIL,INEW,NCALL
      Integer itest,nsum,N12A,N12MX,NDIMO
      LOGICAL FIRST
      EXTERNAL FU,FL,PQCDQQB,PQPM,PQCDI,PQCD,ULMASS

      DATA WDUM/0.D0/,NERP/0/,NX1/10/,NX2/10/,NCMAX/5000/
      DATA DSMALL/1.D-4/,ITEST/1/
      DATA FIRST /.TRUE./
      data nsum/0/
      IF(FIRST) THEN
         FIRST=.FALSE.
         DO K=1,4
            DO I=1,NY
               DO J=1,NQ
                  WMG(K,I,J)=-9999.D0
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      XEPS = 0.005
      IERR = 1
      SLO = 0.0
      SHI = 1.0
      IGENFLO = IGENFL
      NDIMO = NDIM
      IF(IQCDGRID.EQ.0) THEN
         IHERPYSO = IHERPYS
         IF(ITEST.EQ.1) THEN
            IDIRO = IDIR
            IF(IDIRO.EQ.1) THEN
               IHERPYS = 0
            ELSEIF(IDIRO.EQ.0) THEN
               IHERPYS = 1
            ENDIF
            IGENFL = 0
         ENDIF
         if(idir.EQ.1) THEN
            ID = 1
            xpom = 1.d0
         ELSE
            ID= 2
            xpom = DBLE(XFGKI)
            IF(NPOM.EQ.20.OR.NPOM.EQ.21) ID = 3
         endif
c         write(6,*) ' qcdmix IHERPYS =',IHERPYS,IDIR,IDIRO
C select leading order or hihger order process
         IPRO = 12
         RQPM = SNGL(FXN1(X,WDUM))
c         write(6,*) '2nd FXN1 ',RQPM,x
c         call lulist(1)
         IF(ITEST.EQ.1) THEN
            RQPM= SNGL(PQPM(DBLE(YY),Q2))
         ENDIF
         IF(RQPM.LE.0) THEN
            write(6,*) ' THIS MUST NOT HAPPEN '
            write(6,*) ' WMAX = ',WMAX,' RQPM = ',RQPM
            write(6,*) ' IDIR =',IDIR,' NG,NPOM ',NG,NPOM
            write(6,*) ' X ',X
            STOP
         ENDIF
         RQQB=0.0
         RQQBC=0.0
         RQQBB=0.0
         RQQBT=0.0
         RQCDC=0.0
c         write(6,*) 'ihf,ihfla,nflav ',ihf,ihfla,nflav
         IF(IHF.EQ.0) THEN
            IPRO = 13
            AM(1) = 0.0D0
            AM(2) = AM(1)
c calculate maximum of diff x section
C........ gamma gluon fusion
            IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
               WMAX = 0.D0
c            write(6,*) 'qcdmix: IPRO =',IPRO,KPA
               CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQQB)
c            write(6,*) 'qqbar gadap2',RQQB
               WMAX13 = 4.0D0 * WMAX
c           write(6,*) ' gamma gluon fusion'
c           CALL DULIST(1)
            ELSE
cc              write(6,*) xeps,slo,shi
               CALL GADAP(SLO,SHI,PQCDI,XEPS,RQQB)
            ENDIF
c            write(6,*) 'qqbar gadap ',RQQB
            IF(IFPS.NE.10) THEN
C........ QCDC
               IPRO = 15
               IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
c calculate maximum of diff x section
                  WMAX = 0.D0
c           write(6,*) 'qcdmix: IPRO =',IPRO,KPA
                  CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQCDC)
c           write(6,*) 'QCDC gadap2 ',RQCDC
                  WMAX15 = 4.0D0 * WMAX
c           write(6,*) ' QCD compton'
c           CALL DULIST(1)
               ELSE
                  CALL GADAP(SLO,SHI,PQCDI,XEPS,RQCDC)
c           write(6,*) 'QCDC gadap  ',RQCDC
               ENDIF
            ENDIF
         ENDIF
         IF(NFLAV.GE.4.AND.IHF.EQ.0) THEN
C........ gamma gluon fusion charm
            IPRO = 14
c calculate maximum of diff x section
            WMAX = 0.D0
            WMAX14 = 4.0D0 * WMAX
            AM(1) = DBLE(ULMASS(4))
            AM(2) = AM(1)
            IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
c calculate maximum of diff x section
               CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQQBC)
            ELSE
               CALL GADAP(SLO,SHI,PQCDI,XEPS,RQQBC)
            ENDIF
            IF(NFLAV.GE.5) THEN
               AM(1) = DBLE(ULMASS(5))
               AM(2) = AM(1)
               IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
c calculate maximum of diff x section
                  CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQQBB)
               ELSE
                  CALL GADAP(SLO,SHI,PQCDI,XEPS,RQQBB)
               ENDIF
            ENDIF
c            write(6,*) ' rqqbc,rqqbb',rqqbc,rqqbb
         ENDIF
         IF(NFLAV.GE.4.AND.IHF.GT.0) THEN
C........ gamma gluon fusion charm
            IF(IHFLA.EQ.4) THEN
               IPRO = 14
c calculate maximum of diff x section
               WMAX = 0.D0
               WMAX14 = 4.0D0 * WMAX
               AM(1) = DBLE(ULMASS(4))
               AM(2) = AM(1)
               IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                  CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQQBC)
               ELSE
                  CALL GADAP(SLO,SHI,PQCDI,XEPS,RQQBC)
               ENDIF
cc            write(6,*) 'ccbar gadap  RQQBC',RQQBC
            ELSEIF(IHFLA.EQ.5) THEN
               IPRO = 14
               AM(1) = DBLE(ULMASS(5))
               AM(2) = AM(1)
               IF(IQ2.EQ.3.OR.(IQ2.GE.5.AND.IQ2.LE.7)) THEN
                  CALL GADAP2(SLO,SHI,FL,FU,PQCD,XEPS,RQQBB)
               ELSE
                  CALL GADAP(SLO,SHI,PQCDI,XEPS,RQQBB)
               ENDIF
cc            write(6,*) 'bbbar gadap  RQQBB',RQQBB
            ENDIF
         ENDIF
c         write(6,*) ' ihf ',ihf,nflav
c         write(6,*) 'rqpm, rqqb,rqqbc,rqcdc',rqpm,rqqb,rqqbc,rqcdc
c         write(6,*) Q2
         IF(IHF.EQ.0) THEN
c inculde 0.5% uncertainty because of numerics
            RQPMT = RQPM + 0.005*RQPM
            IF(RQPM.LT.(RQQB+RQCDC+RQQBC+RQQBB).AND.
     +       RQPMT.GE.(RQQB+RQCDC+RQQBC+RQQBB)) THEN
               RQPM = RQQB+RQCDC+RQQBC+RQQBB
            ENDIF
            IF(RQPMT.LT.(RQQB+RQCDC+RQQBC+RQQBB)) THEN
               IF(NERP.LT.10) THEN
                  write(6,*) ' ****************************************'
     +            //'*****'
                  write(6,*) ' **  RQQB+RQQBC+RQCDC+RQQBB > RQPM       '
     +            //'   *'
                  write(6,*) ' * RQPM = ',RQPM
                  write(6,*) ' * RQCDC = ',RQCDC,'             *'
                  write(6,*) ' * RQQB = ',RQQB,'               *'
                  write(6,*) ' * RQQBC = ',RQQBC,' RQQBB = ',RQQBB,'*'
                  write(6,*) ' * increase PT2CUT(13):',PT2CUT(13),'  *'
                  write(6,*) ' * increase PT2CUT(15):',PT2CUT(15),'  *'
                  write(6,*) ' * NG =',NG,' NPOM= ',NPOM,' IDIR= ',
     +            IDIR
                  write(6,*) ' ****************************************'
     +            //'*****'
               ELSEIF(NERP.EQ.10) THEN
                  write(6,*) ' ****************************************'
     +            //'*****'
                  write(6,*) ' **         RQQB+RQCDC > RQPM            '
     +            //'    *'
                  write(6,*) ' **         last message printed         '
     +            //'    *'
                  write(6,*) ' ****************************************'
     +            //'*****'
               ENDIF
               NERP = NERP + 1
               RQQB = 0.
               RQQBC = 0.
               RQCDC = 0.
            ENDIF
         ENDIF
c               write(6,*) ' * RQQBC = ',RQQBC,' RQPM = ',RQPM,'*'
         IHERPYS = IHERPYSO
      ELSE
c         write(6,*) 'qcdmix :yy,Q2 ',yy,Q2,ny,nq
c         write(6,*) 'qcdmix : qy ',qy
c         write(6,*) 'qcdmix : qq ',qq
c         write(6,*) ' QCDMIX IDIR = ',IDIR
         IY = 0
   10    IY = IY + 1
         IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
            GOTO 20
         ENDIF
         IF(DBLE(YY).GT.QY(IY+1)) GOTO 10
   20    CONTINUE
         IQ = 0
   30    IQ = IQ + 1
         IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
            GOTO 40
         ENDIF
         IF(Q2.GT.QQ(IQ+1)) GOTO 30
   40    CONTINUE
         IPO = 0
   50    IPO = IPO + 1
         IF(IPO+1.GE.NPPO) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
            GOTO 60
         ENDIF
         IF(DBLE(XFGKI).GT.QPOM(IPO+1)) GOTO 50
   60    CONTINUE
c         write(6,*) 'qcdmix iy,iq',iy,iq,yy,q2
         YD = (DBLE(YY) - QY(IY))/(QY(IY+1) - QY(IY))
         QD = (Q2 - QQ(IQ))/(QQ(IQ+1) - QQ(IQ))
         IF(IDIR.EQ.1) THEN
            X1P=  (QQB(IY+1,IQ)-QQB(IY,IQ))* YD +QQB(IY,IQ)
            X2P=  (QQB(IY+1,IQ+1)-QQB(IY,IQ+1))* YD
     +      +QQB(IY,IQ+1)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            X1P= (QQBDF(IY+1,IQ,IPO)-QQBDF(IY,IQ,IPO))* YD
     +      +QQBDF(IY,IQ,IPO)
            X2P= (QQBDF(IY+1,IQ+1,IPO)-QQBDF(IY,IQ+1,IPO))* YD +
     +      QQBDF(IY,IQ+1,IPO)
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            X1P= (QQBPI(IY+1,IQ,IPO)-QQBPI(IY,IQ,IPO))* YD
     +      +QQBPI(IY,IQ,IPO)
            X2P= (QQBPI(IY+1,IQ+1,IPO)-QQBPI(IY,IQ+1,IPO))* YD +
     +      QQBPI(IY,IQ+1,IPO)
         ENDIF
         RQQB = SNGL((X2P-X1P)*QD + X1P)

         IF(IDIR.EQ.1) THEN
            X1P=  (QQBH(IY+1,IQ)-QQBH(IY,IQ))* YD +QQBH(IY,IQ)
            X2P=  (QQBH(IY+1,IQ+1)-QQBH(IY,IQ+1))* YD +QQBH(IY,
     +      IQ+1)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            X1P= (QQBHDF(IY+1,IQ,IPO)-QQBHDF(IY,IQ,IPO))* YD +
     +      QQBHDF(IY,IQ,IPO)
            X2P= (QQBHDF(IY+1,IQ+1,IPO)-QQBHDF(IY,IQ+1,IPO))* YD +
     +      QQBHDF(IY,IQ+1,IPO)
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            X1P= (QQBHPI(IY+1,IQ,IPO)-QQBHPI(IY,IQ,IPO))* YD +
     +      QQBHPI(IY,IQ,IPO)
            X2P= (QQBHPI(IY+1,IQ+1,IPO)-QQBHPI(IY,IQ+1,IPO))* YD +
     +      QQBHPI(IY,IQ+1,IPO)
         ENDIF
         RQQBC = SNGL((X2P-X1P)*QD + X1P)
         IF(IDIR.EQ.1) THEN
            X1P=  (QQBB(IY+1,IQ)-QQBB(IY,IQ))* YD +QQBB(IY,IQ)
            X2P=  (QQBB(IY+1,IQ+1)-QQBB(IY,IQ+1))* YD +QQBB(IY,
     +      IQ+1)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            X1P= (QQBBDF(IY+1,IQ,IPO)-QQBBDF(IY,IQ,IPO))* YD +
     +      QQBBDF(IY,IQ,IPO)
            X2P= (QQBBDF(IY+1,IQ+1,IPO)-QQBBDF(IY,IQ+1,IPO))* YD +
     +      QQBBDF(IY,IQ+1,IPO)
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            X1P= (QQBBPI(IY+1,IQ,IPO)-QQBBPI(IY,IQ,IPO))* YD +
     +      QQBBPI(IY,IQ,IPO)
            X2P= (QQBBPI(IY+1,IQ+1,IPO)-QQBBPI(IY,IQ+1,IPO))* YD +
     +      QQBBPI(IY,IQ+1,IPO)
         ENDIF
         RQQBB = SNGL((X2P-X1P)*QD + X1P)


         IF(IDIR.EQ.1) THEN
            X1P= (QCDC(IY+1,IQ)-QCDC(IY,IQ))* YD +QCDC(IY,IQ)
            X2P= (QCDC(IY+1,IQ+1)-QCDC(IY,IQ+1))* YD +QCDC(IY,
     +      IQ+1)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            X1P= (QCDCDF(IY+1,IQ,IPO)-QCDCDF(IY,IQ,IPO))* YD +
     +      QCDCDF(IY,IQ,IPO)
            X2P= (QCDCDF(IY+1,IQ+1,IPO)-QCDCDF(IY,IQ+1,IPO))* YD +
     +      QCDCDF(IY,IQ+1,IPO)
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            X1P= (QCDCPI(IY+1,IQ,IPO)-QCDCPI(IY,IQ,IPO))* YD +
     +      QCDCPI(IY,IQ,IPO)
            X2P= (QCDCPI(IY+1,IQ+1,IPO)-QCDCPI(IY,IQ+1,IPO))* YD +
     +      QCDCPI(IY,IQ+1,IPO)
         ENDIF
         RQCDC = SNGL((X2P-X1P)*QD + X1P)
         IF(IFPS.EQ.10) RQCDC  = 0.0
         IF(IDIR.EQ.1) THEN
            X1P= (QPM(IY+1,IQ)-QPM(IY,IQ))* YD +QPM(IY,IQ)
            X2P= (QPM(IY+1,IQ+1)-QPM(IY,IQ+1))* YD +QPM(IY,IQ+
     +      1)
         ELSEIF(IDIR.EQ.0.AND.NPOM.NE.20.AND.NPOM.NE.21) THEN
            X1P= (QPMDF(IY+1,IQ,IPO)-QPMDF(IY,IQ,IPO))* YD
     +      +QPMDF(IY,IQ,IPO)
            X2P= (QPMDF(IY+1,IQ+1,IPO)-QPMDF(IY,IQ+1,IPO))* YD +
     +      QPMDF(IY,IQ+1,IPO)
         ELSEIF(IDIR.EQ.0.AND.(NPOM.EQ.20.OR.NPOM.EQ.21)) THEN
            X1P= (QPMPI(IY+1,IQ,IPO)-QPMPI(IY,IQ,IPO))* YD
     +      +QPMPI(IY,IQ,IPO)
            X2P= (QPMPI(IY+1,IQ+1,IPO)-QPMPI(IY,IQ+1,IPO))* YD +
     +      QPMPI(IY,IQ+1,IPO)
         ENDIF

         RQPM = SNGL((X2P-X1P)*QD + X1P)
c         write(6,*) 'rqpm, rqqb,rqqbc,rqcdc',rqpm,rqqb,rqqbc,rqcdc
         IPRO = 12
         WPA = FXN1(X,WDUM)
         IGENFL = 0
      ENDIF

      IF(IHF.EQ.1) THEN
         IF(IHFLA.EQ.4) RQPM=RQQBC
         IF(IHFLA.EQ.5) RQPM=RQQBB
         IF(IHFLA.EQ.6) RQPM=RQQBT
      ENDIF
c      write(6,*) ' QCDMIX: RQPM,RQQB,RQQBC,RQQBB,RQCDC',
c     + RQPM,RQQB,RQQBC,RQQBB,RQCDC
      PQQB=DBLE(RQQB/RQPM)
      PQQBC=DBLE(RQQBC/RQPM)
      PQQBB=DBLE(RQQBB/RQPM)
      PQCDC=DBLE(RQCDC/RQPM)

      INEW = 0
   70 CONTINUE
c      write(6,*) ' qcdmix pqqb,pqqbc,pqqbb,pqcdc',pqqb,pqqbc,pqqbb,pqcdc
      IF((PQQB+PQCDC+PQQBC+PQQBB).GT.1.02d0) Then
         nsum = nsum + 1
         write(6,*) ' qcdmix: sum >= 1 :IPRO=12 selected'
         write(6,*) ' qcdmix: pqqb,pqqbc,pqqbb ',pqqb,pqqbc,pqqbb
         write(6,*) ' qcdmix: pqcdc ',pqcdc
         write(6,*) ' qcdmix: pqqb+pqqbc+pqqbb+pqcdc ', PQQB+PQCDC+
     +   PQQBC+PQQBB

cc         endif
         IPRO = 12
      else
         PTEST = draprn()
c select process QPM or gamma g fusion
c         write(6,*) 'before sel '
         IF((PQQB+PQCDC+PQQBC+PQQBB).LT.PTEST) THEN
            IPRO = 12
c         write(6,*) ' qpm event '
         ELSEIF((PQQB+PQCDC+PQQBB).LT.PTEST) THEN
            IPRO = 14
            KPA = 4
            WMAX = WMAX14
c         write(6,*) ' cc event '
         ELSEIF((PQQB+PQCDC).LT.PTEST) THEN
            IPRO = 14
            KPA = 5
            WMAX = WMAX14
c         write(6,*) ' bb event '
         ELSEIF(PQQB.LT.PTEST) THEN
            IPRO = 15
            WMAX = WMAX15
c         write(6,*) ' qcdc event '
         ELSE
            IPRO = 13
            WMAX = WMAX13
c         write(6,*) ' qq event '
         ENDIF
c      write(6,*) 'qcdmix y,q2,x_bj',yy,q2,q2/yy/296/296
c      write(6,*) ' ipro ptest ',ipro,ptest
      endif
      IF(IPRO.EQ.12) THEN
         WPA = FXN1(X,WDUM)
      ELSEIF(IPRO.NE.12) THEN
c        IF(IHERAC.EQ.1.and.IDIR.EQ.0) NDIMEN = 5

         IF(IHERAC.EQ.1.and.IDIR.EQ.0) Then
c          NDIM = 3
            NDIMEN = 5
         Endif
         WMAX = 0.D0
         IY = 0
   80    IY = IY + 1
         IF(IY+1.GE.NY) THEN
c          write(6,*) ' y  out of  grid',yy,NY,IY
            GOTO 90
         ENDIF
         IF(DBLE(YY).GT.QY(IY+1)) GOTO 80
   90    CONTINUE
         IQ = 0
  100    IQ = IQ + 1
         IF(IQ+1.GE.NQ) THEN
c          write(6,*) ' q2 out of  grid',q2,NQ,IQ
            GOTO 110
         ENDIF
         IF(Q2.GT.QQ(IQ+1)) GOTO 100
  110    CONTINUE
         IPG = IPRO - 11
c search for maximum of x section
         NNX1 = NX1
         NNX2 = NX2
         IFAIL = 0
         IF(WMG(IPG,IY,IQ).GT.0.D0) THEN
            WMAX=WMG(IPG,IY,IQ)
         ELSE
  120       CONTINUE
            DO 130 I = 0,NNX1
               DO 130 J = 0,NNX2
                  X1 = DFLOAT(I)/DFLOAT(NNX1)
                  Y1 = DFLOAT(J)/DFLOAT(NNX2)
                  IF(X1.EQ.0.D0) X1= DSMALL
                  IF(X1.EQ.1.D0) X1= 1.D0 - DSMALL
                  IF(Y1.EQ.0.D0) Y1= DSMALL
                  IF(Y1.EQ.1.D0) Y1= 1.D0 - DSMALL

                  IF(IDIR.EQ.0) THEN
C generate additional random numbers for order alpha_s QCD processes
C           X(NDIMEN+1) = phi in PHASE
C           X(NDIMEN+2) = cost in PHASE
C           X(NDIMEN+3) = XP2 = E_part/E_proton
C           X(NIDMEN+4) = phi electron
                     X(NDIMEN+1) = draprn()
                     X(NDIMEN+2) = DBLE(X1)
                     X(NDIMEN+3) = DBLE(Y1)
                     X(NDIMEN+4) = draprn()
c              write(6,*) ' qcdmix ',X(NDIMEN+1),X(NDIMEN+2),
c     &           NDIMEN+1,NDIMEN+2
                  ELSEIF(IDIR.EQ.1) THEN
C generate additional random numbers for order alpha_s QCD processes
C           X(NDIMEN+1) = XP2 = E_part/E_proton
C           X(NIDMEN+2) = phi electron
C           X(NDIMEN+3) = phi in PHASE
C           X(NDIMEN+4) = cost in PHASE
                     X(NDIMEN+1) = DBLE(X1)
                     X(NDIMEN+2) = draprn()
                     X(NDIMEN+3) = draprn()
                     X(NDIMEN+4) = DBLE(Y1)
                  ENDIF
                  WPA = FXN1(X,WDUM)
                  IF(WPA.LE.0.0) then
                     WPA = 0.
c                     write(6,*) 'qcdmix:wpa=0,ipro',wpa,ipro,i,j
                  endif
ctest           IF(I.eq.0.or.j.eq.0) write(6,*) ' i,j=0 WPA=',WPA,i,j
                  IF(WPA.GT.WMAX) THEN
                     xm1 = x(ndimen+1)
                     xm2 = x(ndimen+2)
                     xm3 = x(ndimen+3)
                     xm4 = x(ndimen+4)
                     WMAX = WPA
                  ENDIF
c                  IF(XMIN.GT.XMAX) THEN
c                     WRITE(6,*) 'qcdmix: XMIN  >  XMAX  ',XMIN,XMAX,
c     +               IPRO,KPA
c                     SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+
c     +               AM(2)**2))
c                     XMIN=(SMIN+Q2)/(DBLE(Yy)*SSS)
c                     write(6,*) ' 2nd xmin,smin,y,q2 ',xmin,smin,yy,q2
c                     write(6,*) ' 2nd xfgki ',xfgki
c                     IPRO = 12
c                     WPA = FXN1(X,WDUM)
c                     GOTO 140
c                  ENDIF
  130       CONTINUE
            IF(WMAX.LE.0.D0) THEN
c            WRITE(6,*) ' WMAX always 0.0: old IPRO:',IPRO
               NNX1 = NNX1 + 25
               NNX2 = NNX2 + 25
c               WRITE(6,*) ' increase nr points in  grid ',NNX1,NNX2
               IFAIL = IFAIL + 1
               IF(IFAIL.LE.1) THEN
                  GOTO 120
               ELSE
c                  WRITE(6,*) ' WMAX always 0.0: old IPRO:',IPRO
c                  WRITE(6,*) ' too many trials: select new process '
c                  write(6,*) ' IDISDIF ',IDISDIF,' IDIR ',IDIR
c                  write(6,*) ' IHERPYS ',IHERPYS
                  INEW = INEW + 1
                  IF(INEW.LE.2) THEN
                     GOTO 70
                  ELSE
                     IPRO = 12
                     WPA = FXN1(X,WDUM)
                     write(6,*) ' too many trials: select IPRO =12'
                     write(6,*) ' IDIR,XFGKI,IPO ',IDIR,XFGKI,IPO
                     GOTO 160
                  ENDIF
               ENDIF
            ENDIF
ctest            WMAX = 3.D0 * WMAX
            WMAX = 1.5D0 * WMAX
ccccc            WMG(IPG,IY,IQ)=WMAX
         ENDIF
         NCALL = 0
  140    CALL draprnV(XG,4)
         DO 150 I=1,4
  150    X(NDIMEN+I)=XG(I)
c                     X(NDIMEN+1) = xm1
c                     X(NDIMEN+2) = xm2
c                     X(NDIMEN+3) = xm3
c                     X(NDIMEN+4) = xm4
c              write(6,*) ' qcdmix 2nd',X(NDIMEN+1),X(NDIMEN+2),
c     &           NDIMEN+1,NDIMEN+2
         NCALL = NCALL+1
         WPA = FXN1(X,WDUM)
         WTEST=draprn()
cnew
         IF(WPA.LT.0.0) THEN
            IPRO = 12
            WPA = FXN1(X,WDUM)
            LST(21) = 31
            GOTO 160
         ENDIF
cnew
         IF(WMAX.LE.0.0) RETURN
         IF(WPA.GT.WMAX) THEN
            write(6,*) ' WPA > WMAX for process IPRO = ',IPRO, WPA,
     +      WMAX
            write(6,*) ' in bin IY,IQ,IPG  ',IY,IQ,IPG
c            write(6,*) ' for IDIR = ',IDIR
c            write(6,*) ' X(NDIMEN+I) ', (X(NDIMEN+I),I=1,4)
c            write(6,*) ' max at  ',xm1,xm2,xm3,xm4
c            write(6,*) ' PQQB,PQQBC,PQCDC ',PQQB,PQQBC,PQCDC,PTEST
            WMG(IPG,IY,IQ)=WPA
            IPRO = 12
            WPA = FXN1(X,WDUM)
            LST(21) = 30
         ELSE
c            write(6,*) ' qcdmix: ipro,wpa,wmax ',ipro,wpa,wmax
            IF(WPA.LE.0.0D0.AND.NCALL.LT.NCMAX) THEN

               GOTO 140
            ELSEIF(WPA.LE.0.0D0) THEN
               write(6,*) ' WPA LE 0: switch to IPRO=12',IPRO,WPA
               IPRO = 12
               WPA = FXN1(X,WDUM)
               LST(21) = 31
               GOTO 160
            ENDIF
            IF(WPA/WMAX.LT.WTEST.AND.NCALL.LT.NCMAX) THEN
               GOTO 140
            ELSEIF(NCALL.GE.NCMAX) THEN
               write(6,*) ' NCALL GT 1000: switch to IPRO=12',IPRO
               IPRO = 12
               WPA = FXN1(X,WDUM)
               LST(21) = 32
               GOTO 160
            ENDIF
         ENDIF
      ENDIF
  160 CONTINUE
      IF(WPA.NE.0.0D0) THEN
         IERR = 0
      else
         wtest = FXN1(X,WDUM)
         write(6,*) 'QCDMIX WPA = 0 IPRO ',WPA,x,IPRO,wtest
      ENDIF
      NDIM = NDIMO
      IGENFL = IGENFLO
      RETURN
      END
