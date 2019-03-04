*CMZ :  2.08/04 22/12/99  15.39.26  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  17.08.08  by  Hannes Jung
*CMZ :  2.07/03 23/05/99  14.51.10  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE PRODIFF
      IMPLICIT None
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


*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      DOUBLE PRECISION PARL
      REAL CUT,XLEP,YLEP,W2LEP,Q2LEP,ULEP
	Integer LST
      COMMON /RAPTOU/CUT(14),LST(40),PARL(30),XLEP,YLEP,W2LEP,Q2LEP,ULEP
      Double Precision XMX
      COMMON/TEST/XMX
      DOUBLE PRECISION PHV(4)
      Double Precision draprn
	Real CHI,PTSPL,PHISPL
      REAL SNGL
	Integer I,ipom,kpart,npdc1,npdc2,iprot,KFLCH,KFLSP,KSIGN,KC2
	integer LUCOMP
      common/pdcyn/npdc1,npdc2
	Double Precision R,AMSP,AMK2,beta,A,shr,phr,pe,pz,pw1
	Double Precision pms,pms2,pms4,xmas,fc,dot
      EXTERNAL draprn
      IPROT = -9999
	PTSPL = 0.
	PHISPL = 0.
      DO 10 I=3,N
         IF(IABS(K(I,2)).EQ.2212) IPROT = I
         IF(IABS(K(I,2)).EQ.100) IPOM = I
   10 CONTINUE
      IF(IPROT.LT.0) RETURN
      K(IPROT,2)=2210
      KPART = K(IPROT,2)
      K(IPROT,1)=11
      KSIGN=ISIGN(1,KPART)
      KFLCH=0
C...GLUON REMOVED.
      R=6.D0*draprn()
      IF(R.LT.3.) THEN
         KFLCH=2*KSIGN
         KFLSP=2101*KSIGN
      ELSEIF(R.LT.4.) THEN
         KFLCH=2*KSIGN
         KFLSP=2103*KSIGN
      ELSE
         KFLCH=1*KSIGN
         KFLSP=2203*KSIGN
      ENDIF
      K(N+1,1)=1
      K(N+2,1)=1
      K(N+1,2)= KFLCH
      K(N+2,2)= KFLSP
      K(N+1,3)=IPROT
      K(N+2,3)=IPROT
   20 CONTINUE
C...Split energy-momentum of target remnant according to functions P(z)
C...z=E-pz fraction for qq (q) forming jet-system with struck Q (Qbar)
C...1-z=E-pz fraction for qQbar (qqQ) hadron
C...mq=mass of (light) parton remnant q (qq) in jet system
C...mQ=mass of produced (heavy flavour) hadron
      AMSP=DBLE(ULMASS(KFLSP))
      AMK2=DBLE(ULMASS(KFLCH))
      IF(IREM.EQ.1) THEN
C...Old Lepto treatment
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=1/3 since a=1 fixed
         BETA=DBLE(1.-SQRT(draprn()))
C...Flip if baryon produced
         KC2=IABS(LUCOMP(KFLCH))
         IF(KC2.GE.301.AND.KC2.LE.400) BETA=1.D0-BETA
      ELSEIF(IREM.EQ.2) THEN
C...Update of LEPTO treatment
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=mq/(mq+mQ) --> a=a(mq,mQ)
         A=(AMSP+AMK2)/AMSP - 2.D0
         BETA=draprn()**(1.D0/(A+1.D0))
      ELSEIF(IREM.EQ.3) THEN
C...Using Peterson fragmentation function
C...P(z)=N/(z(1-1/z-c/(1-z))**2)  where c=(mq/mQ)**2  (FC=-c)
         FC=-(AMSP/AMK2)**2
   30    BETA=draprn()
         IF(-4.D0*FC*BETA*(1.D0-BETA)**2.LT.draprn()*
     +    ((1.D0-BETA)**2-FC*BETA)**2) GOTO 30
      ELSEIF(IREM.LT.0) THEN
         beta = DBLE(XFGKI + (SNGL(XF) - XFGKI)*draprn())

      ENDIF
      IF(beta.lt.1.5*xfgki) goto 20
C --
      chi = XFGKI/SNGL(beta)
c use constituent masses in ULMASS
cccc      MSTJ(93)=1
      P(N+1,5)=DBLE(ULMASS(K(N+1,2)))
      P(N+1,1)=P(2,1)*DBLE(CHI) - P(IPOM,1) + DBLE(PTSPL*COS(PHISPL))
      P(N+1,2)=P(2,2)*DBLE(CHI) - P(IPOM,2) + DBLE(PTSPL*SIN(PHISPL))
      P(N+1,1) = P(2,1)*DBLE(CHI) - P(IPOM,1)
      P(N+1,2) = P(2,2)*DBLE(CHI) - P(IPOM,2)
c use constituent masses in ULMASS
cccc      MSTJ(93)=1
      P(N+2,5)=DBLE(ULMASS(K(N+2,2)))

C--
      DO 40 I=1,4
   40 PHV(I) = P(IPOM,I) + P(1,I)
      SHR = SSS
      SHR = SQRT(SHR)
      PHR = DOT(PHV,PHV)
      P(N+2,1) = P(IPROT,1) - P(N+1,1)
      P(N+2,2) = P(IPROT,2) - P(N+1,2)
      PMS2 = P(N+1,5)**2 + P(N+1,1)**2 + P(N+1,2)**2
      PMS4 = P(N+2,5)**2 + P(N+2,1)**2 + P(N+2,2)**2
      PMS = PMS2/DBLE(CHI) +PMS4/(1.D0-DBLE(CHI))
      PE = 0.5D0*(SHR + (PMS - PHR)/SHR)
      PZ = SQRT(PE**2 - PMS)
      PW1=DBLE(CHI)*(PE + PZ)
      P(N+1,4) = 0.5D0*(PW1 + PMS2/PW1)
      P(N+1,3) = 0.5D0*(PW1 - PMS2/PW1)

      P(N+2,4) = PE - P(N+1,4)
      P(N+2,3) = PZ - P(N+1,3)
ccc      xm1 = P(N+1,4)**2 - P(N+1,1)**2-P(N+1,2)**2-P(N+1,3)**2
ccc      xm2 = P(N+2,4)**2 - P(N+2,1)**2-P(N+2,2)**2-P(N+2,3)**2
ccc      xm1=dsqrt(xm1)
ccc      xm2=dsqrt(xm2)
ccc      write(6,*) ' masses xm1 ',xm1,p(n+1,5)
ccc      write(6,*) ' masses xm2 ',xm2,p(n+2,5)

      xmx = DOT1(N+1,N+1) + 2.D0*DOT1(N+1,N+2) + DOT1(N+2,N+2)
      xmas = dsqrt(xmx)
      if(xmas.lt.1.1) then
         npdc1 = 0
         npdc2 = 0
         DO I = 1,5
            P(N+1,I) = P(IPROT,I)
         ENDDO
         K(N+1,1)=1
         K(N+1,2)=2212
         K(N+1,3)=IPROT
         N=N+1
      ELSE
         npdc1 = N+1
         npdc2 = N+2
         K(npdc1,1) = 12
         K(npdc2,1) = 11

         N = N+2

      ENDIF

      RETURN
      END
