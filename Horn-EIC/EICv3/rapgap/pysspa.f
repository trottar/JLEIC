*CMZ :  2.08/04 22/12/99  15.39.26  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  18.07.54  by  Hannes Jung
*CMZ :  2.07/03 23/05/99  18.05.05  by  Hannes Jung
*CMZ :          23/05/99  15.18.18  by  Hannes Jung
C from PYTHIA 4.8 manual
C IPY(14) (D=2) details of initial state parton showers
C         =0    not included
C         >=1   included
C NEW
C IORDER  =0    no ordering
C         =1    Q2 values at branches are strictly ordered, increasing towards
C               the hard scattering
C         =2    Q2 and opening angles of emitted (on shell or time like) partons
C               are both strictly ordered, increasing towards the hard
C               interaction.
C
C IALPS   =0
C         =1    alphas first order with scale Q2
C         =2    alphas first order with scale k_t**2=(1-z)*Q2
C
C ITIMSHR =0    no shower of time like partons
C         =1    time like partons may shower
C ISOFTG (D=1) treatment of soft gluons
C         =0    soft gluons are entirely neglected
C         =1    soft gluons are resummed and included together with the hard
C               radiation as an effective z shift.
C
C PYPAR(23) (D=2) effective minimum energy (in CM frame)  of timelike or
C               on shell parton emitted in spacelike showers.
C PYPAR(24) (D=0.001) effective lower cutoff in 1-z in spacelike showers
C               in addition to PYPAR(23)
*-- Author :    Hannes Jung   03/04/94


C.. THIS IS A COPY FROM LEPTO61
C.  CHANGE for photoproduction and gamma glu fusion

      SUBROUTINE PYSSPA(IPU1,IPU2)
      IMPLICIT None
C...NEW X REDEFINITION
C...GENERATES SPACELIKE PARTON SHOWERS
      REAL CUT,XLP,YLP,W2LP,Q2LP,ULP
      DOUBLE PRECISION PARL
	Integer LLST
      COMMON /RAPTOU/ CUT(14),LLST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
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

*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      REAL PYPAR,PYVAR
	Integer IPY
      COMMON /PYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL X,SH,TH,UH,Q2
	Integer ISUB,KFL
      COMMON /PYPROC/ ISUB,KFL(3,2),X(2),SH,TH,UH,Q2
      REAL XQ
      COMMON /MYINT1/ XQ(2,-6:6)
      REAL XA,XFA,XB,XFB,Q2REF,XBMIN
      DOUBLE PRECISION ROBO,THE
	Integer IFLS
	Double Precision xs,zs,q2s,tevs,xfs,wtap,wtsf
      DIMENSION IFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVS(2),ROBO(5),
     +XFS(2,-6:6),XFA(-6:6),XFB(-6:6),WTAP(-6:6),WTSF(-6:6)
c     +XFS(2,-6:6),XFA(-25:25),XFB(-25:25),WTAP(-6:6),WTSF(-6:6)
      DOUBLE PRECISION DQ2(3),DSH,DSHZ,DSHR,DPLCM,DPC(3),DPD(4),DMS,
     +DMSMA,DPT2,DPB(4),DBE1(4),DBE2(4),DBEP,DGABEP,DPQ(4),DPQS(2),
     +DM2,DQ2B,DROBO(5),DBEZ,THE2(2)
C-GI &DQ23,DPH(4),DM2,DQ2B,DQM2
c.hju
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

*KEND.
      Double Precision draprn
	Real QMAX
      EXTERNAL draprn
	Integer LST,IRES
      COMMON/EPPARA/ LST(30),IRES(2)
      REAL SNGL
	Integer I,J,NS,IFLA,IFLB,NQ,IS,IPU1,IPU2,iremp,isea,ilep,ipo
	Integer i1,i2,itemp,ikin,kn1,kd1,ir,jb,id1,id2
	Integer ntry,ntry2,nprint,ncall,ifl,jt,IHFC,IHFX,IHFT,jr,it
	Double Precision z,xe,xe0,q2e,tmax,b0,q2b,tevb,qmass,wtapq,wtsum
	Double Precision VALENCE,SEARN,WTRAN,WTZ,XBNEW,RSOFT,ZU,ALPRAT
	Double Precision THE2T
      LOGICAL LAST
      LOGICAL LPRINT,FIRST
      LOGICAL OLDEPS
*KEEP,PSHWR.
      INTEGER IORDER,IALPS,ITIMSHR,ISOFTG,ICCFM
      COMMON /PSHWR/IORDER,IALPS,ITIMSHR,ISOFTG,ICCFM
*KEND.
      DATA LPRINT/.TRUE./
      DATA FIRST/.TRUE./
      DATA OLDEPS/.FALSE./
      DATA NPRINT /0/
      DATA NCALL /0/
      DATA IFLA,NQ/0,0/,Z,XE0,XA/3*0./,DSHZ,DMSMA,DPT2,DSHR/4*0.D0/
      DATA IS/2*0/
c set flag for epsilon treatment
c      OLDEPS = .TRUE.
cccccccccccccccccccccccccccccccc
      NTRY = 0
      NCALL = NCALL + 1
      IF(NCALL.GE.NPRINT) LPRINT=.FALSE.
c set flag for new remnant treatment (motivated by LEPTO 6.5)
      IREMP = 1
      LAST=.FALSE.
      ISEA = 0
c no special sea quark treatment for pomeron
      IF(IPY(42).EQ.100) IREMP = 0
c no special sea quark treatment for resolved gammas
      IF(IPY(41).EQ.22) IREMP = 0
C...COMMON CONSTANTS, SET UP INITIAL VALUES
      ILEP=0
      IF(IPU1.EQ.0) ILEP=1
      IF(IPU2.EQ.0) ILEP=2
      IF(ILEP.eq.1.and.ILEPTO.eq.0) ILEP=0
cccc
      PYPAR(23)=0.0
check pypar
cccc      PYPAR(23)=4.0
cccc
      IF(ILEP.EQ.0) PYPAR(23)=2.0
      DO 10  I=1,3
   10 DQ2(I) = 0.D0
      SH=SNGL(P(25,5)**2)
      IF(N.GE.27) SH=SNGL(P(27,5)**2)
      IF(IPRO.EQ.99) SH=SNGL(P(28,5)**2)
      CALL LSCALE(-1,QMAX)
      Q2E=DBLE(QMAX**2)
      TMAX=DLOG(Q2E/DBLE(PYPAR(21)**2))
      IF(LPRINT) THEN
         write(6,*) ' pysspa TMAX,PYPAR(26),PYPAR(27),Q2E,PYPAR(21)',
     +   TMAX,PYPAR(26),PYPAR(27),Q2E,PYPAR(21)
      ENDIF


      IF(ILEP.GE.1) THEN

         SH=SNGL(P(25,5)**2)
         IF(N.GE.27) SH=SNGL(P(27,5)**2)
         IF(IPRO.EQ.99) SH=SNGL(P(28,5)**2)
         CALL LSCALE(-1,QMAX)
         Q2E=DBLE(QMAX**2)
c.hju        Q2E=MAX(PYPAR(21)**2,MIN(Q2E,(0.95/X(3-ILEP)-1.)*Q2-SH,
c.hju     &  Q2/2.+SH))
c.hju        Q2E=MAX(PYPAR(21)**2,MIN(Q2E,(0.95/X(3-ILEP)-1.)*Q2LP-SH,
c.hju     &  Q2LP/2.+SH))
c.hju        Q2E=MAX(PYPAR(21)**2,MIN(Q2E,(1.-x(3-ilep))*PARL(22),
c.hju     &  Q2LP/2.+SH))
         Q2E=DMAX1(DBLE(PYPAR(21)**2), DMIN1(Q2E,DBLE(-Q2LP+PYVAR(31)*
     +   SNGL(PARL(21))-SH), DBLE(Q2LP/2.+SH)))
         TMAX=DLOG(Q2E/DBLE(PYPAR(21))**2)
      ENDIF
      IF(LPRINT) THEN
         write(6,*) ' start pysspa ',PYPAR(26)*sngl(Q2E),
     +   MAX(PYPAR(22),2.*PYPAR(21)**2)
         write(6,*) 'Q2E ',Q2E,' SH ',SH,' Q2LP ',Q2LP,PYVAR(31),
     +   PARL(21)
      ENDIF
      IF(SNGL(Q2E).LT.MAX(PYPAR(22),2.*PYPAR(21)**2).OR.
     +TMAX.LT.0.2) THEN
         IF(LPRINT) THEN
            write(6,*) ' PYYSPA 1st check: '
            write(6,*) ' PYPAR(26)=',PYPAR(26),' PYPAR(22)',PYPAR(22)
            write(6,*) ' PYPAR(21) ',PYPAR(21)
            write(6,*) ' SH ',SH,' Q2LP ',Q2LP,PYVAR(31),PARL(21)
            write(6,*) ' -Q2LP+PYVAR(31)*PARL(21)-SH', -Q2LP+PYVAR(31)*
     +      sngl(PARL(21))-SH
            write(6,*) ' Q2E ',Q2E,' TMAX ',TMAX
            write(6,*) ' QMAX ',QMAX,' SH ',SH,IPRO
         ENDIF
         LST(21)=55
         RETURN
      ENDIF
      IF(LPRINT) write(6,*) 'pysspa:check ok '
      IF(ILEP.EQ.0) XE0=DBLE(2.*PYPAR(23)/PYVAR(1))
      B0=DBLE((33.-2.*IPY(8))/6.)
      NS=N
      MSTU(2)=0
   20 N=NS
      IF(ILEP.GE.1) THEN
         NQ=IPU2-2
         IF(ILEP.EQ.2) NQ=IPU1+2
         DPQS(1)=P(NQ,3)
         DPQS(2)=P(NQ,4)
         IF(IPRO.EQ.12) THEN
            XBMIN=X(3-ILEP)*MAX(0.5,SH/Q2LP)
            CALL PYSTFU(IPY(43-ILEP),XBMIN,Q2,XFB)
cccccc        write(6,*) ' XFB = ',(XFB(i),i=-6,6)
            DO 30  IFL=-6,6
   30       XQ(3-ILEP,IFL)=XFB(IFL)
         ENDIF
      ENDIF
      DO 40  JT=1,2
         IFLS(JT)=KFL(2,JT)
         IF(KFL(2,JT).EQ.21) IFLS(JT)=0
         IFLS(JT+2)=IFLS(JT)
         XS(JT)=DBLE(X(JT))
         ZS(JT)=1.D0
         IF(ILEP.EQ.0) Q2S(JT)=Q2E
         THE2(JT)=100.D0
         TEVS(JT)=TMAX
ccc         write(6,*) ' TEVS,JT,TMAX',TEVS(JT),JT,TMAX
         DO 40 IFL=-6,6
c.. checking x values HJU
cccc            IF(XS(JT).GT.0.8.and.Ilep.ne.0) XQ(JT,IFL) = 0.
            IF(XS(JT).GT.0.8) XQ(JT,IFL) = 0.
c switch off ps from photon....
cccc            IF(ILEP.EQ.0) XQ(1,IFL) = 0.
c.. end checking
   40 XFS(JT,IFL)=DBLE(XQ(JT,IFL))
      IF(ILEP.GE.1) THEN
         Q2S(ILEP)=P(NQ,5)**2
         DQ2(ILEP)=Q2S(ILEP)
c         write(6,*) 'Q2E ',Q2E
         Q2S(3-ILEP)=Q2E
      ENDIF
cadded
      IF(ISEMIH.EQ.1) THEN
         Q2S(2)=P(IPU2,5)**2
         DQ2(2)=Q2S(2)
ccc       write(6,*) ' pysspa ISEMIH:',DQ2(2)
      ENDIF
      DSH=DBLE(SH)
      IHFC=0
      IHFX=0

C...PICK UP LEG WITH HIGHEST VIRTUALITY
   50 CONTINUE
      IF(N.GT.MSTU(4)-10) THEN
         WRITE(6,*) ' PYSSPA: no more memory in LUJETS'
         LST(21)=51
         RETURN
      ENDIF
      DO 60  I=N+1,N+8
         DO 60 J=1,5
            K(I,J)=0
   60 P(I,J)=0.D0
C     CALL GULIST(21,2)
      N=N+2
      JT=1
      IF((N.GT.NS+2.AND.Q2S(2).GT.Q2S(1).AND.ILEP.EQ.0).OR.ILEP.EQ.1)
     +JT=2
      IF(LPRINT) THEN
         write(6,*) ' pysspa JT,Q2S(1),Q2S(2) ',JT,Q2S(1),Q2S(2)
      ENDIF
      JR=3-JT
      IFLB=IFLS(JT)
      XB=SNGL(XS(JT))
c.hju      IF(ILEP.GE.1.AND.N.EQ.NS+2) XB=XS(JT)*MAX(SH/Q2,0.5)
      IF(ILEP.GE.1.AND.N.EQ.NS+2.AND.IPRO.EQ.12)
     +                     XB=SNGL(XS(JT))*MAX(SH/Q2LP,0.5)
      DO 70  IFL=-6,6
   70 XFB(IFL)=SNGL(XFS(JT,IFL))

      Q2B=Q2S(JT)
      TEVB=TEVS(JT)
      IF(IORDER.EQ.0) THEN
         IF(ZS(JT).GT.0.9999) THEN
            Q2B=Q2S(JT)
         ELSE
            Q2B=0.5D0*(1.D0/ZS(JT)+1.D0)*Q2S(JT)+0.5D0*(1.D0/ZS(JT)-
     +      1.D0) *(Q2S(3-JT)- DSH+DSQRT((DSH+Q2S(1)+Q2S(2))**2+8.D0*
     +      Q2S(1)*Q2S(2)* ZS(JT)/(1.D0-ZS(JT))))
         ENDIF
         TEVB=DLOG(Q2B/DBLE(PYPAR(21)**2))
      ENDIF
      IF(ILEP.EQ.0) THEN
         DSHR=2.D0*DSQRT(DSH)
         DSHZ=DSH/DBLE(ZS(JT))
      ELSEIF(ILEP.GE.1) THEN
         DSHZ=DSH
         IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/ZS(JT)-DQ2(JR)+
     +   DBLE(PYPAR(22))
         DPD(2)=DSHZ+DQ2(JR)+DBLE(PYPAR(22))

         QMASS=DBLE(ULMASS(IABS(IFLB)))
         IF(IABS(IFLB).EQ.0) QMASS=DBLE(ULMASS(21))
C...CHECK IF QUARK PAIR CREATION ONLY POSSIBILITY
         IF(DQ2(JR).LE.4.D0*QMASS**2) THEN
            DM2=QMASS**2
            DPC(1)=DQ2(JR)*(DBLE(PYPAR(22))+DM2)**2
            DPC(2)=DPD(2)*(DPD(2)-2D0*DBLE(PYPAR(22)))*
     +       (DBLE(PYPAR(22))+DM2)
            DPC(3)=DBLE(PYPAR(22))*(DPD(2)-2D0*DBLE(PYPAR(22)))**2
c.hju here are few changes for Q2 = 0
            IF(Q2LP.GE.0.5*PYPAR(22)) THEN
               XE0=1D0-(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/
     +         (2D0*DPC(1))
            ELSE
               XE0=1.D0-(DBLE(PYPAR(22))*(DPD(2)-2D0*DBLE(PYPAR(22)))
     +         /DPD(2)/(DBLE(PYPAR(22))+DM2))
            ENDIF
            if(oldeps) then
            else
               XE0 = dble(xb)*xe0/(1.d0 - xe0)
            endif
c.hju end of changes
         ELSE
c.hju here are few changes for Q2 = 0
            IF(Q2LP.GE.0.5*PYPAR(22)) THEN
               XE0=1D0-(DPD(2)-2D0*DBLE(PYPAR(22)))*(DPD(2)-
     +         DSQRT(DPD(2)**2- 4D0*DQ2(JR)*DBLE(PYPAR(22))))/(2D0*
     +         DQ2(JR)*DBLE(PYPAR(22)))
            ELSE
               XE0=1.D0-(DBLE(PYPAR(22))*(DPD(2)-2D0*DBLE(PYPAR(22)))
     +          /DPD(2)/DBLE(PYPAR(22)))
            ENDIF
            if(oldeps) then
            else
               XE0 = dble(xb)*xe0/(1.d0 - xe0)
            endif
         ENDIF
         if(oldeps) then
         else
c include fixed cut off to avoid soft gluons
c also solves energy momentum conservation problem
            XE0=DMAX1(XE0,2.D0*DBLE(PYPAR(23))/DSQRT(PARL(21)))
         endif
      ENDIF
   80 XE=DMAX1(XE0,DBLE(XB)*(1.d0/(1.d0-DBLE(PYPAR(24)))-1.d0))
      IF(XB+SNGL(XE).GE.0.999) THEN
         Q2B=0.D0
         IF(LPRINT) THEN
            write(6,*) ' here goto 150 ',XB,XE
         ENDIF
         GOTO 180
      ENDIF
      IF(LPRINT) THEN
         write(6,*) 'calc A P weights ',JT
      ENDIF
C...CALCULATE ALTARELLI-PARISI AND STRUCTURE FUNCTION WEIGHTS
      DO 90  IFL=-6,6
         WTAP(IFL)=0.D0
   90 WTSF(IFL)=0.D0
      IF(IFLB.EQ.0) THEN
         WTAPQ=DBLE(16.*(1.-SQRT(XB+SNGL(XE)))/(3.*SQRT(XB)))
         DO 100 IFL=-IPY(8),IPY(8)
            IF(IFL.EQ.0) WTAP(IFL)=6.D0*DLOG((1.D0-DBLE(XB))/XE)
  100    IF(IFL.NE.0) WTAP(IFL)=WTAPQ
      ELSE
         WTAP(0)=0.5D0*DBLE(XB)*(1.D0/(DBLE(XB)+XE)-1.D0)
         WTAP(IFLB)=8.D0*DLOG((1.D0-DBLE(XB))*(DBLE(XB)+XE)/XE)/3.D0
      ENDIF
  110 NTRY = NTRY+1
      IF(NTRY.GT.500) THEN
         LST(21)=56
         RETURN
      ENDIF
      WTSUM=0.D0
      IF(IHFC.EQ.0) THEN
         DO 120 IFL=-IPY(8),IPY(8)
            IF(LPRINT) THEN
               write(6,*) ' in loop IFL,IFLB,XFB(IFL)',IFL,IFLB,
     +         XFB(IFL)
            ENDIF
            WTSF(IFL)=DBLE(XFB(IFL)/MAX(1E-10,XFB(IFLB)))
  120    WTSUM=WTSUM+WTAP(IFL)*WTSF(IFL)
         IF(IABS(IFLB).GE.4.AND.WTSUM.GT.1E3) THEN
            IHFX=1
            DO 130 IFL=-IPY(8),IPY(8)
  130       WTSF(IFL)=WTSF(IFL)*1D3/WTSUM
            WTSUM=1D3
         ENDIF
      ENDIF

C...CHOOSE NEW T AND FLAVOUR
      IF(LPRINT) THEN
         write(6,*) ' before : TEVB,B0,WTSUM',TEVB,B0,WTSUM
      ENDIF
      NTRY2 = 0
  140 NTRY2 = NTRY2 + 1
      IF(NTRY2.GT.500) THEN
         LST(21)=56
         RETURN
      ENDIF
      IF(IALPS.EQ.1) THEN
         TEVB =TEVB*draprn()**(B0/DMAX1(0.0001d0,WTSUM))
      ELSEIF(IALPS.EQ.2) THEN
         TEVB =TEVB*draprn()**(B0/DMAX1(0.0001d0,5.D0*WTSUM))
      ELSE
         write(6,*) ' fatal IALPS not selected ',IALPS
         write(6,*) ' PROGRAM STOPPED '
         STOP
      ENDIF
c translate t into Q2 scale
      Q2REF=PYPAR(21)**2*EXP(SNGL(TEVB))
      IF(LPRINT) THEN
         write(6,*) 'TEVB ',TEVB,' WTSUM ',WTSUM,' Q2REF ',Q2REF
      ENDIF
      Q2B=DBLE(Q2REF)
      DQ2B=Q2B
      IF(ILEP.GE.1) THEN
         DSHZ=DSH
         IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/DBLE(ZS(JT))-DQ2(JR)+
     +   DQ2B
      ENDIF
      ISEA = 0
      IF(LPRINT) THEN
         write(6,*) ' pysspa IPY(42)',IPY(42)
      ENDIF
      IF(IPY(40+JT).EQ.22) THEN
check that Q2B > Q2LP (virtuality of photon) for resolved gammas
         IF(Q2B.LE.Q2LP) THEN
            IF(LPRINT) THEN
               write(6,*) 'res.gamma  Q2B,Q2LP ',Q2B,Q2LP
            ENDIF
            Q2B=0.D0
         ENDIF
      ENDIF
changed      IF(Q2B.LT.PYPAR(22).AND.IREMP.EQ.1.AND.LAST) THEN
      IF(Q2B.LT.PYPAR(22).AND.IREMP.EQ.1.AND.LAST.AND.
     +IPY(40+JT).NE.22) THEN
         IREMP = 0
c check if last parton was valence or sea parton
         IF(IFLB.NE.0) THEN
            VALENCE =  DBLE(XFB(IFLB)-XFB(-IFLB))
            IF(LPRINT) THEN
               write(6,*) ' PYSSPA: IFLB,XFB(IFLB),VALENCE', IFLB,
     +         XFB(IFLB),XFB(-IFLB),VALENCE
            ENDIF
            SEARN = draprn()
            IF(LPRINT) THEN
               write(6,*) ' PYSSPA SEA ',SEARN,sngl(VALENCE)/XFB(IFLB)
            ENDIF
            IF(VALENCE.LT.1.D-5) THEN
               ISEA = 1
            ELSE
               IF(SEARN.GT.DABS(VALENCE/DBLE(XFB(IFLB)))) ISEA = 1
            ENDIF
         ENDIF
         IF(ISEA.EQ.1) THEN
            IF(LPRINT) THEN
               write(6,*) ' pysspa ISEA = ',ISEA
            ENDIF
            IFLA = 0
            Q2B = DBLE(PYPAR(22))
            DQ2B = Q2B
            IF(ILEP.GE.1) THEN
               DSHZ=DSH
               IF(N.GT.NS+4) DSHZ=(DSH+DQ2(JR)-DQ2(JT))/DBLE(ZS(JT))-
     +         DQ2(JR)+DQ2B
            ENDIF
         ENDIF
      ENDIF
      IF(Q2B.LT.PYPAR(22).AND.ISEA.EQ.0) THEN
         IF(LPRINT) THEN
            write(6,*) ' Q2B =',Q2B,PYPAR(22),ISEA
         ENDIF
         Q2B=0.D0
      ELSE
         IF(ISEA.EQ.1) GOTO 160
         LAST=.TRUE.
         IF(LPRINT) THEN
            write(6,*) '2nd. JT = ',JT,Q2B
         ENDIF
         WTRAN=draprn()*WTSUM
         IFLA=-IPY(8)-1
  150    IFLA=IFLA+1
         WTRAN=WTRAN-WTAP(IFLA)*WTSF(IFLA)
         IF(IFLA.LT.IPY(8).AND.WTRAN.GT.0.) GOTO 150
  160    CONTINUE
C...CHOOSE Z VALUE AND CORRECTIVE WEIGHT
c... g --> g + g.
         IF(IFLB.EQ.0.AND.IFLA.EQ.0) THEN
            Z=DBLE(1./(1.+((1.-XB)/XB)*(SNGL(XE)/(1.-XB))**draprn()))
            WTZ=(1.D0-Z*(1.D0-Z))**2
c... q --> g + q
         ELSEIF(IFLB.EQ.0) THEN
            Z=DBLE(XB/(1.-draprn()*(1.-SQRT(XB+SNGL(XE))))**2)
            WTZ=0.5D0*(1.D0+(1.D0-Z)**2)*DSQRT(Z)
c... g --> q + q_bar
         ELSEIF(IFLA.EQ.0) THEN
            Z=DBLE(XB*(1.+draprn()*(1./(XB+SNGL(XE))-1.)))
            WTZ=1.D0-2.D0*Z*(1.D0-Z)
c... q --> q + g.
         ELSE
            Z=DBLE(1.-(1.-XB)*(SNGL(XE)/((XB+SNGL(XE))*
     +       (1.-XB)))**draprn())
            WTZ=0.5D0*(1.D0+Z**2)
         ENDIF
c         pause 'after z '
C...REWEIGHT FIRST LEG BECAUSE OF MODIFIED XB OR CHECK PHASE SPACE
         IF(ILEP.GE.1.AND.N.EQ.NS+2) THEN
            IF(IPRO.EQ.12) THEN
               XBNEW=DBLE(X(JT))*(1.D0+(DSH-Q2B)/DQ2(JR))
            ELSE
c. new calculation for first x_g rescaling
               XBNEW=DBLE(X(JT))*(1.D0 + (DQ2(JR)-Q2B)/DSH)
            ENDIF
            IF(XBNEW.GT.MIN(Z,0.999d0)) GOTO 140
            XB=SNGL(XBNEW)
         ENDIF
         IF(ISEA.EQ.1) GOTO 170
C...SUM UP SOFT GLUON EMISSION AS EFFECTIVE Z SHIFT
         IF(ISOFTG.GE.1) THEN
            RSOFT=6.D0
            IF(IFLB.NE.0) RSOFT=8.D0/3.D0
            Z=Z*(TEVB/TEVS(JT))**(RSOFT*XE/((DBLE(XB)+XE)*B0))
            IF(Z.LE.XB) GOTO 140
         ENDIF
C...CHECK IF HEAVY FLAVOUR BELOW THRESHOLD
         IHFT=0
         IF(ILEP.GE.1.AND.IABS(IFLB).GE.4.AND.(XFB(IFLB).LT.1E-10.OR.
     +   Q2B.LT.5.*ULMASS(IABS(IFLB))**2)) THEN
            IHFT=1
            IFLA=0
         ENDIF
c include this also for res gammas
         IF(ILEP.EQ.0.AND.IABS(IFLB).GE.4.AND.(XFB(IFLB).LT.1E-10.OR.
     +   Q2B.LT.5.*ULMASS(IABS(IFLB))**2)) THEN
            IHFT=1
            IFLA=0
         ENDIF

C...FOR LEPTOPRODUCTION, CHECK Z AGAINST NEW LIMIT
         IF(ILEP.GE.1) THEN
            DPD(2)=DSHZ+DQ2(JR)+DQ2B
            DM2=DBLE(ULMASS(IABS(IFLA-IFLB))**2)
            IF(IABS(IFLA-IFLB).EQ.0) DM2=DBLE(ULMASS(21)**2)
            DPC(1)=DQ2(JR)*(DQ2B+DM2)**2
            DPC(2)=DPD(2)*(DPD(2)-2D0*DQ2B)*(DQ2B+DM2)
            DPC(3)=DQ2B*(DPD(2)-2D0*DQ2B)**2
c.hju here are few changes for Q2 = 0
            IF(Q2LP.GE.0.5*PYPAR(22)) THEN
               ZU=(DPC(2)-DSQRT(DPC(2)**2-4D0*DPC(1)*DPC(3)))/(2D0*
     +         DPC(1))
            ELSE
               ZU = DQ2B*(DPD(2)-2D0*DQ2B)/DPD(2)/(DQ2B+DM2)
            ENDIF
c.end of changes for Q2 = 0
c            write(6,*) ' Z limit ',Z,ZU
            IF(Z.GE.ZU) GOTO 140
         ENDIF
changed
C...OPTION WITH EVOLUTION IN KT2=(1-Z)Q2:
changed         IF(IPY(14).GE.5.AND.IPY(14).LE.6.AND.N.LE.NS+4) THEN
C...CHECK THAT (Q2)LAST BRANCHING < (Q2)HARD
changed            IF(LPRINT) write(6,*) ' check branching '
changed            IF(Q2B/(1.D0-Z).GT.PYPAR(26)*Q2) GOTO 140
changed         ELSEIF(IPY(14).GE.3.AND.IPY(14).LE.6.AND.N.GE.NS+6) THEN
C...CHECK THAT Z,Q2 COMBINATION IS KINEMATICALLY ALLOWED
changed            IF(LPRINT) write(6,*) 'check z,Q2 combination'
changed            Q2MAX=0.5D0*(1.D0/ZS(JT)+1.D0)*DQ2(JT)+0.5D0*
changed     +       (1.D0/ZS(JT)-1.D0)*
changed     +      (DQ2(3-JT)-DSH+SQRT((DSH+DQ2(1)+DQ2(2))**2+
changed     +      8.D0*DQ2(1)*DQ2(2)
changed     +      * ZS(JT)/(1.D0-ZS(JT))))
changed            IF(Q2B/(1.D0-Z).GE.Q2MAX) GOTO 140
changed         ENDIF
         IF(IALPS.GE.2) THEN
C...OPTION WITH ALPHAS((1-Z)Q2): DEMAND KT2 > CUTOFF, REWEIGHT
            IF((1.D0-Z)*Q2B.LT.PYPAR(22)) GOTO 140
            ALPRAT=TEVB/(TEVB+DLOG(1.D0-Z))
            IF(ALPRAT.LT.5.*draprn()) GOTO 140
            IF(ALPRAT.GT.5.) WTZ=WTZ*ALPRAT/5.D0
         ENDIF
check for angular ordering
         IF(IORDER.EQ.2.AND.NTRY2.LT.200) THEN
            THE2T=(4.D0*Z**2*Q2B)/(dble(PYVAR(1))*(1.D0-Z)*dble(XB)**2)
            IF(THE2T.GT.THE2(JT)) GOTO 140
         ENDIF
C...WEIGHTING WITH NEW STRUCTURE FUNCTIONS
         IF(LPRINT) THEN
            write(6,*) 'CALL PYSTFU(IPY(40+JT),XB,Q2REF,XFB)',IPY(40+
     +      JT)
            write(6,*) ' Q2REF,Q2LP ',Q2REF,Q2LP
         ENDIF
         CALL PYSTFU(IPY(40+JT),XB,Q2REF,XFB)
         XA=XB/SNGL(Z)
         IF(LPRINT) THEN
            write(6,*) 'CALL PYSTFU(IPY(40+JT),XA,Q2REF,XFA)',IPY(40+
     +      JT)
         ENDIF
         CALL PYSTFU(IPY(40+JT),XA,Q2REF,XFA)
         IF(IHFT.EQ.1.OR.IHFX.EQ.1) THEN
            IF(XFA(IFLA).LT.1E-10) IHFC=1
            IF(LPRINT) THEN
               write(6,*) 'IFLA = ',IFLA,' xa ',xa,ipy(40+JT), XFA(IFLA
     +         ),Q2REF
               write(6,*) ' goto 180 '
            ENDIF
            GOTO 180
         ELSEIF(XFB(IFLB).EQ.0.0) THEN
            IF(LPRINT) THEN
               write(6,*) 'IFLB = ',IFLB,' xb ',xb,ipy(40+JT), XFB(IFLB
     +         ),Q2REF
               write(6,*) ' return '
            ENDIF
            LST(21) = 54
c            pause
            RETURN
         ELSEIF(XFB(IFLB).LT.1E-20) THEN
            IF(LPRINT) THEN
               write(6,*) ' goto 110 '
            ENDIF
            GOTO 110
         ENDIF
         IF(LPRINT) THEN
            write(6,*) 'WTZ,XFA(IFLA)/XFB(IFLB),WTSF(IFLA)', WTZ,
     +      XFA(IFLA)/XFB(IFLB),WTSF(IFLA),IFLA,IFLB
         ENDIF
         IF(SNGL(WTZ)*XFA(IFLA)/XFB(IFLB).LT.
     +    draprn()*SNGL(WTSF(IFLA))) THEN
            IF(LPRINT) THEN
               write(6,*) 'WTZ check ',N,NS+2
               write(6,*) ' goto 80 or 110'
            ENDIF
            IF(ILEP.GE.1.AND.N.EQ.NS+2) GOTO 80
            IF(LPRINT) THEN
               write(6,*) 'goto 110'
            ENDIF
            GOTO 110
         ENDIF
c         pause 'after loop'
  170    CONTINUE
      ENDIF
      IF(LPRINT) THEN
         write(6,*) 'before HARD scatter in cm frame ',N,NS+4-2*MIN(1,
     +   ILEP)
      ENDIF

  180 IF(N.EQ.NS+4-2*MIN(1,ILEP)) THEN
         IF(LPRINT) THEN
            write(6,*) 'HARD scatter in cm frame '
         ENDIF
C...DEFINE TWO HARD SCATTERERS IN THEIR CM-FRAME
         DQ2(JT)=Q2B
changed         IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.D0-Z)
         IF(ILEP.EQ.0) THEN
            DPLCM=DSQRT((DSH+DQ2(1)+DQ2(2))**2-4.D0*DQ2(1)*DQ2(2))/DSHR
            DO 190 JR=1,2
               I=NS+2*JR-1
               IPO=19+2*JR
               K(I,1)=14
               K(I,2)=IFLS(JR+2)
               IF(IFLS(JR+2).EQ.0) K(I,2)=21
               IF(IFLS(JR+2).EQ.22) K(I,1)=11
               K(I,3)=0
               K(I,4)=IPO
               K(I,5)=IPO
               IF(K(IPO,2).EQ.22) K(I,4)=0
               IF(K(IPO,2).EQ.22) K(I,5)=0
               P(I,1)=0.D0
               P(I,2)=0.D0
               P(I,3)=DPLCM*(-1)**(JR+1)
               P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR
               P(I,5)=-SQRT(DQ2(JR))
               K(I+1,1)=-1
               K(I+1,2)=K(IPO+1,2)
               K(I+1,3)=I
               K(I+1,4)=0
               K(I+1,5)=0
               P(I+1,1)=0.D0
               P(I+1,2)=0.D0
               P(I+1,3)=IPO
               P(I+1,4)=IPO
               P(I+1,5)=0.D0
               P(IPO+1,1)=I
               P(IPO+1,2)=I
               K(IPO,4)=MOD(K(IPO,4),MSTU(5))+I*MSTU(5)
               K(IPO,5)=MOD(K(IPO,5),MSTU(5))+I*MSTU(5)
               IF(K(IPO,2).EQ.22) K(I,4)=0
               IF(K(IPO,2).EQ.22) K(I,5)=0
  190       CONTINUE
         ELSE
C..LEPTOPRODUCTION EVENTS: BOSON AND HADRON REST FRAME
            I1=NS+2*ILEP-1
            I2=NS-2*ILEP+5
            DO 200 ITEMP=NS+1,NS+4
               DO 200 J=1,5
                  K(ITEMP,J)=0
  200       P(ITEMP,J)=0.D0
            DO 210 J=1,5
  210       P(I1,J)=P(NQ,J)
            K(NS+1,1)=11
            K(NS+3,1)=14
            IF(ILEP.EQ.2) THEN
               K(NS+1,1)=14
               K(NS+3,1)=11
            ENDIF
            K(NS+2,1)=-1
            K(NS+4,1)=-1
            K(NS+1,3)=0
            K(NS+2,3)=NS+1
            K(NS+3,3)=0
            K(NS+4,3)=NS+3
            K(I1,2)=KFL(2,ILEP)
            K(I2,2)=KFL(2,3-ILEP)
            DPD(1)=DSH+DQ2(1)+DQ2(2)
            DPD(3)=(3-2*ILEP)*DSQRT(DPD(1)**2-4D0*DQ2(1)*DQ2(2))
c.hju here are few changes for Q2 = 0
            IF(Q2LP.GE.0.5*PYPAR(22)) THEN
               P(I2,3)=(DPQS(2)*DPD(3)-DPQS(1)*DPD(1))/ (2D0*DQ2(JR))
               P(I2,4)=(DPQS(1)*DPD(3)-DPQS(2)*DPD(1))/ (2D0*DQ2(JR))
            ELSE
               P(I2,3) =- (DPD(1)**2 + 4.d0*DPQS(2)**2 *DQ2(2))/
     +         (4.d0*DPQS(1)*DPD(1))
               P(I2,4) = (DPD(1)**2 - 4.d0*DPQS(1)**2 *DQ2(2))/ (4.d0*
     +         DPQS(2)*DPD(1))
               P(I2,3)=-DPD(1)/(2D0*(DPQS(1)+DPQS(2)))-
     +                 DQ2(2)*DPQS(2)/DPD(1)
               P(I2,4)=DPD(1)/(2D0*(DPQS(1)+DPQS(2)))-
     +                 DQ2(2)*DPQS(1)/DPD(1)
cc               write(6,*) ' pysspa new ',P(I2,3),P(23,3)
cc               write(6,*) ' pysspa new',P(I2,4),P(23,4)
Cend test
               IF(P(I2,4).LT.0.D0) THEN
c                  write(6,*) 'P(I2,4) lt 0.:'
c                  goto 20
               endif
            ENDIF
c.hju end of these changes
            P(I2,5)=-SQRT(DQ2(3-ILEP))
            P(I2+1,3)=MAX(IPU1,IPU2)
            P(I2+1,4)=MAX(IPU1,IPU2)
            K(I2,4)=K(I2,4)-MOD(K(I2,4),MSTU(5))+MAX(IPU1,IPU2)
            K(I2,5)=K(I2,5)-MOD(K(I2,5),MSTU(5))+MAX(IPU1,IPU2)
            P(26-2*ILEP,1)=I2
            P(26-2*ILEP,2)=I2
            K(25-2*ILEP,4)=MOD(K(25-2*ILEP,4),MSTU(5))+I2*MSTU(5)
            K(25-2*ILEP,5)=MOD(K(25-2*ILEP,5),MSTU(5))+I2*MSTU(5)
            N=N+2
         ENDIF

      ELSEIF(N.GT.NS+4) THEN
         IF(LPRINT) THEN
            write(6,*) 'find max allowed mass'
         ENDIF
C...FIND MAXIMUM ALLOWED MASS OF TIMELIKE PARTON
changed         IF(ILEP.EQ.0) JR=3-JT
         DQ2(3)=Q2B
changed         IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(3)=Q2B/(1.D0-Z)
         IF(IS(1).GE.1.AND.IS(1).LE.MSTU(4)) THEN
            DPC(1)=P(IS(1),4)
            DPC(3)=0.5D0*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))
         ELSE
C...IS(1) not initialized
            DPC(1)=0.D0
            DPC(3)=0.5D0*( 0.D0 +ABS(P(IS(2),3)))
         ENDIF
         DPC(2)=P(IS(2),4)
         DPD(1)=DSH+DQ2(JR)+DQ2(JT)
         DPD(2)=DSHZ+DQ2(JR)+DQ2(3)
         DPD(3)=DSQRT(DPD(1)**2-4.D0*DQ2(JR)*DQ2(JT))
         DPD(4)=DSQRT(DPD(2)**2-4.D0*DQ2(JR)*DQ2(3))
         IKIN=0
c.hju        IF((Q2S(JR).GE.0.5*PYPAR(22).AND.DPD(1)-DPD(3).GE.1D-10*DPD(1))
c.hju     &  .OR.ILEP.GE.1) IKIN=1
         IF((Q2S(JR).GE.0.5*PYPAR(22).AND.DPD(1)-DPD(3).GE.1D-10*DPD(1)
     +   )) IKIN=1
         IF(IKIN.EQ.0.AND.ILEP.EQ.0) DMSMA=(DQ2(JT)/DBLE(ZS(JT))-DQ2(3)
     +   ) *(DSH/(DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))
         IF(IKIN.EQ.0.AND.ILEP.GE.1) DMSMA=(DPD(1)**2 * DQ2(3) +
     +   DPD(2)**2 * DQ2(JT))/DPD(1)/DPD(2) -DQ2(JT)-DQ2(3)
         IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/(2.D0*DQ2(JR)
     +   )-DQ2(JT)-DQ2(3)
         IF(DMSMA.LT.0.0) THEN
            IF(LPRINT) THEN
               WRITE(6,*) 'DMSMA = ',DMSMA,' IKIN = ',IKIN,' JT = ',JT
               write(6,*) 'DPC(2)=',DPC(2),' DPD(1)=',DPD(1)
               write(6,*) 'DPD(2)=',DPD(2),' DPD(3)=',DPD(3)
               write(6,*) 'DPD(4)=',DPD(4),' DQ2(JT)=',DQ2(JT),JT
               write(6,*) 'DQ2(JR)=',DQ2(JR),JR
               write(6,*) 'ZS(JT)=',ZS(JT),' DQ2(3)=',DQ2(3)
               write(6,*) 'DSH= ',DSH,'DSHZ=',DSHZ
               write(6,*) ' XS(JT) = ',XS(JT),X(JT)
            ENDIF
         ENDIF
C...GENERATE TIMELIKE PARTON SHOWER (IF REQUIRED)
         IT=N-1
         K(IT,1)=3
         K(IT,2)=IFLB-IFLS(JT+2)
         IF(IFLB-IFLS(JT+2).EQ.0) K(IT,2)=21
         P(IT,5)=DBLE(ULMASS(K(IT,2)))
         IF(DMSMA.LE.P(IT,5)**2) THEN
            IF(LPRINT) THEN
               write(6,*) 'DMSMA LT P(IT,5)**2',DMSMA,P(IT,5)**2
               write(6,*) ' really goto 20 ????? '
            ENDIF
            GOTO 20
         ENDIF
         P(IT,2)=0.D0
         DO 220 J=1,5
            K(IT+1,J)=0
  220    P(IT+1,J)=0.D0
         K(IT+1,1)=-1
         K(IT+1,2)=K(IS(JT)+1,2)
         K(IT+1,3)=IT
         IF(ITIMSHR.EQ.1) THEN
            P(IT,1)=0.D0
            IF(ILEP.EQ.0) P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR
            IF(ILEP.GE.1) P(IT,4)=0.5D0*(P(IS(JT),3)*DPD(2)+ DPQS(1)*
     +      (DQ2(JT)+DQ2(3)+P(IT,5)**2))/(P(IS(JT),3)*DPQS(2)- P(IS(JT)
     +      ,4)*DPQS(1))-DPC(JT)
c            IF(P(IT,4).LT.0.0) write(6,*) 'P(IT,4) < 0.:'
            P(IT,3)=DSQRT(DMAX1(0.D0,P(IT,4)**2-P(IT,5)**2))
            CALL DUSHOW(IT,0,SQRT(MIN(SNGL(DMSMA),PYPAR(25)*Q2)))
            IF(N.GE.IT+2) P(IT,5)=P(IT+2,5)
            IF(N.GT.MSTU(4)-10) THEN
               WRITE(6,*) ' PYSSPA: no more memory in LUJETS'
               LST(21)=52
               RETURN
            ENDIF
            DO 230 I=N+1,N+8
               DO 230 J=1,5
                  K(I,J)=0
  230       P(I,J)=0.D0
         ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: TIMELIKE PARTON SHOWER
         DMS=P(IT,5)**2
c.hju
         DPT2 = -1.D0
c.hju
         IF(IKIN.EQ.0.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/
     +   (DSH+DQ2(JT))
         IF(IKIN.EQ.1.AND.ILEP.EQ.0) DPT2=(DMSMA-DMS)*(0.5D0*DPD(1)*
     +   DPD(2)+0.5D0*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     +   (4.D0*DSH*DPC(3)**2)
c.test        IF(IKIN.EQ.1.AND.ILEP.GE.1) DPT2=(DMSMA-DMS)*(0.5*DPD(1)*
         IF(ILEP.GE.1) DPT2=(DMSMA-DMS)*(0.5D0*DPD(1)* DPD(2)+
     +    0.5D0*DPD(3)*
     +   DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/ DPD(3)**2
c.hju
c.test        IF(IKIN.EQ.0.AND.ILEP.GE.1) DPT2=(DMSMA-DMS)*DPD(2)/DPD(1)
c.hju
         IF(DPT2.LT.0.) THEN
            WRITE(6,*) 'DPT2 ',DPT2
            write(6,*) ' really goto 20 ????? '
            GOTO 20
         ENDIF
         K(IT,3)=N+1
         P(IT,1)=SQRT(DPT2)
         IF(ILEP.EQ.0) THEN
            DPB(1)=(0.5D0*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/
     +      DSHR)/DPC(3)-DPC(3)
            P(IT,3)=DPB(1)*(-1)**(JT+1)
            P(IT,4)=(DSHZ-DSH-DMS)/DSHR
         ELSE
            DPC(3)=DQ2(JT)+DQ2(3)+DMS
            DPB(2)=DPQS(2)*P(IS(JT),3)-DPQS(1)*DPC(JT)
            DPB(1)=0.5D0*(DPC(JT)*DPD(2)+DPQS(2)*DPC(3))/DPB(2)-
     +      P(IS(JT),3)
            P(IT,3)=DPB(1)
            P(IT,4)=0.5D0*(P(IS(JT),3)*DPD(2)+ DPQS(1)*DPC(3))/DPB(2)-
     +      DPC(JT)
ctest for numeric stability
            P(IT,4) = DSQRT(P(IT,3)**2+P(IT,1)**2+P(IT,5)**2)
cend test
            IF(P(IT,4).LE.0.0) THEN
               write(6,*) 'pysspa P(4)<0 ',
     +       0.5D0*(P(IS(JT),3)*DPD(2)+ DPQS(1)*DPC(3))/DPB(2),DPC(JT)
c               pause
c               CALL DULIST(1)
            ENDIF
         ENDIF
         IF(N.GE.IT+2) THEN
            MSTU(1)=IT+2
            DPB(1)=DSQRT(DPB(1)**2+DPT2)
            DPB(2)=DSQRT(DPB(1)**2+DMS)
            DPB(3)=P(IT+2,3)
            DPB(4)=DSQRT(DPB(3)**2+DMS)
            DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)*
     +      DPB(1))
            CALL DUDBRB(MSTU(1),MSTU(2),0.D0,0.D0,0.D0,0.D0,DBEZ)
            THE=DLANGL(P(IT,3),P(IT,1))
            CALL DUDBRB(MSTU(1),MSTU(2),THE,0.D0,0.D0,0.D0,0.D0)
c test
            MSTU(1) = 0
cend test
         ENDIF

C...RECONSTRUCT KINEMATICS OF BRANCHING: SPACELIKE PARTON
         K(N+1,1)=14
         K(N+1,2)=IFLB
         IF(IFLB.EQ.0) K(N+1,2)=21
         K(N+1,3)=0
         P(N+1,1)=P(IT,1)
         P(N+1,2)=0.D0
         P(N+1,3)=P(IT,3)+P(IS(JT),3)
         P(N+1,4)=P(IT,4)+P(IS(JT),4)
         P(N+1,5)=-SQRT(DQ2(3))

c check that p_t of cascade is smaller than p_t of hard process
         IF(P(N+1,1).GE.SQRT(PT2H)) THEN
c           write(6,*) ' warning from pysspa: '
c        write(6,*) ' p_t of branching > p_t of hard scattering '
c     +  ,P(N+1,1),SQRT(PT2H),P(N+1,5)
         ENDIF
         DO 240 J=1,5
            K(N+2,J)=0
  240    P(N+2,J)=0.D0
         K(N+2,1)=-1
         K(N+2,2)=K(IS(JT)+1,2)
         K(N+2,3)=N+1

C...DEFINE COLOUR FLOW OF BRANCHING
         K(IS(JT),1)=14
         K(IS(JT),3)=N+1
         ID1=IT
         KN1=ISIGN(500+IABS(K(N+1,2)),2*K(N+1,2)+1)
         KD1=ISIGN(500+IABS(K(ID1,2)),2*K(ID1,2)+1)
         IF(K(N+1,2).EQ.21) KN1=500
         IF(K(ID1,2).EQ.21) KD1=500
         IF((KN1.GE.501.AND.KD1.GE.501).OR.(KN1.LT.0.AND. KD1.EQ.500)
     +   .OR.(KN1.EQ.500.AND.KD1.EQ.500.AND. draprn().GT.0.5)
     +   .OR.(KN1.EQ.500.AND.KD1.LT.0)) ID1=IS(JT)
         ID2=IT+IS(JT)-ID1
         P(N+2,3)=ID1
         P(N+2,4)=ID2
         P(ID1+1,1)=N+1
         P(ID1+1,2)=ID2
         P(ID2+1,1)=ID1
         P(ID2+1,2)=N+1
         K(N+1,4)=K(N+1,4)-MOD(K(N+1,4),MSTU(5))+ID1
         K(N+1,5)=K(N+1,5)-MOD(K(N+1,5),MSTU(5))+ID2
         K(ID1,4)=MOD(K(ID1,4),MSTU(5))+(N+1)*MSTU(5)
         K(ID1,5)=MOD(K(ID1,5),MSTU(5))+ID2*MSTU(5)
         K(ID2,4)=MOD(K(ID2,4),MSTU(5))+ID1*MSTU(5)
         K(ID2,5)=MOD(K(ID2,5),MSTU(5))+(N+1)*MSTU(5)
         N=N+2
C...BOOST TO NEW CM-FRAME
         MSTU(1)=NS+1
         IF(ILEP.EQ.0) THEN
            CALL DUDBRB(MSTU(1),MSTU(2),0.D0,0.D0,
     +      -(P(N-1,1)+P(IS(JR),1))/(P(N-1,4)+P(IS(JR),4)),
     +      0.D0,-(P(N-1,3)+P(IS(JR),3))/(P(N-1,4)+P(IS(JR),4)))
            IR=N-1+(JT-1)*(IS(1)-N+1)
            CALL DUDBRB(MSTU(1),MSTU(2), -DLANGL(P(IR,3),P(IR,1)),
     +      DBLE(PARU(2)*draprn()),0.D0,0.D0,0.D0)
         ELSE
C...REORIENTATE EVENT WITHOUT CHANGING THE BOSON FOUR MOMENTUM
            DO 250 J=1,4
  250       DPQ(J)=P(NQ,J)
            DBE1(4)=DPQ(4)+P(N-1,4)
            DO 260 J=1,3,2
  260       DBE1(J)=-(DPQ(J)+P(N-1,J))/DBE1(4)
c          write(6,*) ' DBE1(1)=',DBE1(1),'DBE1(3)=',DBE1(3)
            IF((1D0-DBE1(1)**2-DBE1(3)**2).LE.0.0) THEN
               MSTU(1) = 0
               LST(21) = 53
               RETURN
            ENDIF
            DBE1(4)=1D0/DSQRT(1D0-DBE1(1)**2-DBE1(3)**2)
            DBEP=DBE1(1)*DPQ(1)+DBE1(3)*DPQ(3)
            DGABEP=DBE1(4)*(DBE1(4)*DBEP/(1D0+DBE1(4))+DPQ(4))

            DO 270 J=1,3,2
  270       DPQ(J)=DPQ(J)+DGABEP*DBE1(J)
            DPQ(4)=DBE1(4)*(DPQ(4)+DBEP)
            DPC(1)=DSQRT(DPQ(1)**2+DPQ(3)**2)
            DBE2(4)=-(DPQ(4)*DPC(1)-DPQS(2)*DSQRT(DPQS(2)**2+DPC(1)**2-
     +      DPQ(4)**2))/(DPC(1)**2+DPQS(2)**2)
            THE=DLANGL(DPQ(3),DPQ(1))
            DBE2(1)=DBE2(4)*SIN(THE)
            DBE2(3)=DBE2(4)*COS(THE)
            DBE2(4)=1D0/(1D0-DBE2(1)**2-DBE2(3)*2)

C...CONSTRUCT THE COMBINED BOOST
            DPB(1)=DBE1(4)**2*DBE2(4)/(1D0+DBE1(4))
            DPB(2)=DBE1(1)*DBE2(1)+DBE1(3)*DBE2(3)
            DPB(3)=DBE1(4)*DBE2(4)*(1D0+DPB(2))
            DO 280 JB=1,3,2
  280       DROBO(JB+2)=(DBE1(4)*DBE2(4)*DBE1(JB)+DBE2(4)*DBE2(JB)+
     +      DPB(1)*DBE1(JB)*DPB(2))/DPB(3)
            IF(DSQRT(DROBO(3)**2+DROBO(5)**2).GT.0.99999999D0) THEN
               MSTU(1) = 0
               LST(21) = 53
               RETURN
            ENDIF
            CALL DUDBRB(MSTU(1),MSTU(2),0.D0,0.D0,
     +      DROBO(3),0.D0,DROBO(5))
            IF(ILEP.EQ.1) THE=DLANGL(P(NS+1,3),P(NS+1,1))
            IF(ILEP.EQ.2) THE=DBLE(PARU(1))+DLANGL(P(NS+3,3),P(NS+3,1))
            CALL DUDBRB(MSTU(1),MSTU(2),-THE, DBLE(PARU(2)*draprn()),
     +      0D0,0D0,0D0)
         ENDIF
         MSTU(1)=0
      ENDIF

C...SAVE QUANTITIES, LOOP BACK
      IF(LPRINT) THEN
         write(6,*) ' save quantities and loop back '
      ENDIF
      IS(JT)=N-1
      IF(ILEP.EQ.2.AND.N.EQ.NS+4) IS(JT)=N-3
      Q2S(JT)=Q2B
      DQ2(JT)=Q2B
      IF(IORDER.EQ.2) THE2(JT)=THE2T
changed      IF(IPY(14).GE.3.AND.IPY(14).LE.6) DQ2(JT)=Q2B/(1.D0-Z)
      DSH=DSHZ
      IF(Q2B.GE.0.5*PYPAR(22)) THEN
         IFLS(JT+2)=IFLS(JT)
         IFLS(JT)=IFLA
         XS(JT)=DBLE(XA)
         ZS(JT)=Z
         DO 290 IFL=-6,6
  290    XFS(JT,IFL)=DBLE(XFA(IFL))
         TEVS(JT)=TEVB
      ELSE
         IF(JT.EQ.1) IPU1=N-1
         IF(JT.EQ.2) IPU2=N-1
      ENDIF
      IF(LPRINT) THEN
         write(6,*) 'check end: N=',N,' NS+2=',NS+2
         write(6,*) MAX(IABS(1-ILEP)*Q2S(1),MIN(1,2-ILEP)*Q2S(2)),
     +   0.5*PYPAR(22)
         write(6,*) Q2S(1),Q2S(2)
      ENDIF
      IF(MAX(IABS(1-ILEP)*Q2S(1),MIN(1,2-ILEP)*Q2S(2)).GE.0.5*PYPAR(22)
     +.OR.N.LE.NS+2) GOTO 50
      IF(ILEP.EQ.0) THEN
C...BOOST HARD SCATTERING PARTONS TO FRAME OF SHOWER INITIATORS
         DO 300 J=1,3
  300    DROBO(J+2)=(P(NS+1,J)+P(NS+3,J))/(P(NS+1,4)+P(NS+3,4))
         DO 310 J=1,5
  310    P(N+2,J)=P(NS+1,J)
         MSTU(1)=N+2
         MSTU(2)=N+2
         CALL DUDBRB(N+2,N+2,0.D0,0.D0,-DROBO(3),-DROBO(4),-DROBO(5))
         ROBO(2)=DLANGL(P(N+2,1),P(N+2,2))
         ROBO(1)=DLANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
c.test        MSTU(1)=4
         MSTU(1)=19
         MSTU(2)=NS
         CALL DUDBRB(MSTU(1),MSTU(2), ROBO(1),ROBO(2),DROBO(3),DROBO(4)
     +   ,DROBO(5))
         MSTU(1)=0
         MSTU(2)=0
      ENDIF

C...STORE USER INFORMATION
      K(21,1)=14
c.test      IF(ILEP.NE.0) K(21,1)=11
      IF(IRES(1).EQ.0) K(21,1)=11
      K(23,1)=14
      K(21,3)=NS+1
      K(23,3)=NS+3
      DO 320 JT=1,2
         KFL(1,JT)=IFLS(JT)
         IF(IFLS(JT).EQ.0) KFL(1,JT)=21

  320 PYVAR(30+JT)=SNGL(XS(JT))
      DO 330 I=NS+1,N
         DO 330 J=1,5
  330 V(I,J)=0.
c      CALL DULIST(1)
      RETURN
      END
