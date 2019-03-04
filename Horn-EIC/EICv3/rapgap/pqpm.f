*CMZ :  2.08/02 14/07/99  18.40.55  by  Hannes Jung
*CMZ :  2.08/01 23/06/99  07.56.53  by  Hannes Jung
*CMZ :  2.08/00 17/06/99  17.21.39  by  Hannes Jung
*CMZ :  2.07/03 23/05/99  18.42.36  by  Hannes Jung
*-- Author :
      FUNCTION PQPM(YY1,Q2)
      IMPLICIT None
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


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

*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEND.
      REAL XPQ(-6:6)
      REAL SNGL,XS,Q2S
      REAL SCALE
      Integer KPART
      COMMON/PINT/ SCALE,KPART
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
      Double Precision YY1,Q2,GEV2NB,X,F2,ALPH_EM,SIGMA,WTGLU,GF,WMAT
      Double Precision PQPM,beta
      Integer I,LUCHGE
      DATA GEV2NB/.3893857D+6/

      X = Q2/YY1/SSS
      XS = REAL(X)
      Q2S = REAL(Q2)
c      CALL PYSTFU(K(2,2),SNGL(X),SNGL(Q2),XPQ)
      do i=-6,6
         XPQ(i)=0.0
      enddo
      IF(ID.EQ.1) THEN
         CALL PYSTFU(2212,SNGL(X),SNGL(Q2),XPQ)
      ELSE
         BETA = X/XPOM
c     write(6,*) xpom,x,beta,id
         IF(BETA.LT.1.and.beta.gt.0.) THEN
c            write(6,*) xpom,x,beta,id
            CALL RASTFU(KINT(2,2),SNGL(beta),SNGL(Q2),XPQ)
c        write(6,*) KINT(2,2),SNGL(beta),SNGL(Q2),XPQ
         ENDIF
      ENDIF
c check for exclusive heavy flavor production
      if(IHF.GT.0) then
         do i=-6,6
            if(iabs(i).ne.IHFLA) xpq(i)=0
         enddo
      endif
c      write(6,*) ' PQPM x,q2,xpq     ',ihf,ihfla,x,q2,(xpq(i),i=0,5)

ccc      write(6,*) ' PQPM x,q2,xpq MXcut',x,q2,(xpq(i),i=0,4)
      F2 = 0.D0
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
      IF(INTER.LT.2) THEN
         DO 10 I=-NFLAV,NFLAV
   10    F2 = F2 + DBLE(XPQ(I))*DFLOAT(LUCHGE(I))**2/9.D0
cccc         F2 = 2*DBLE(XPQ(4))*DFLOAT(LUCHGE(4))**2/9.D0
         SIGMA = 2.D0*PI*ALPH_EM**2 /YY1 /Q2/Q2 *(1.d0+(1.d0-yy1)**2)*
     +   F2
      ELSEIF(INTER.EQ.2) THEN
         WTGLU = 0.D0
         DO 20  I=1,NFLAV-1,2
            IF(ISIGN(1,K(1,2)).EQ.1) THEN
               WTGLU = WTGLU + (1.D0 - YY1 + YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ( -I)+XPQ(I+1)) + (YY1 - YY1**2/2.D0)*2.D0*
     +         DBLE(-XPQ(-I)+ XPQ(I+1))
            ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
               WTGLU = WTGLU + (1.D0 - YY1 + YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ(I )+XPQ(-I-1)) - (YY1 - YY1**2/2.D0)*2.D0*
     +         DBLE(XPQ(I)-XPQ(- I-1))
            ENDIF
   20    CONTINUE
         GF = PI*ALPH_EM/(SIN2W*XMW2*DSQRT(2.D0))
         WMAT = GF**2/2.D0/PI /(1.D0 + Q2/XMW2)**2 /YY1
capply additional factor 0.5 because only right/or left handed electrons contr.
         WMAT = WMAT * 0.5D0
         SIGMA = WMAT * WTGLU
c         write(6,*) ' here in CC pqpm',sigma*gev2nb
      ENDIF
c      write(6,*) pi,alph,q2,f2,yy1,q2
      PQPM = SIGMA * GEV2NB
      RETURN
      END
