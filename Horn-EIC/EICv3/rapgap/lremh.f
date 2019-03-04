*CMZ :  2.08/04 22/12/99  15.39.25  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  21.54.41  by  Hannes Jung
*-- Author :
C **********************************************************************

      SUBROUTINE LREMH(IFLRO,IFLR,K2,Z)
	Implicit None

C...Gives flavour and energy-momentum fraction Z for the particle
C...to be produced out of the target remnant when that is not a
C...simple diquark.
      DOUBLE PRECISION PARL
	Real CUT,X,Y,W2,Q2,U
	Integer LST
      COMMON /RAPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      Double Precision draprn
	Real ULMASS
	Integer IFLRO,IFLQ,IFLQQ,IFLR,KSP,K2,IDUM,K2A,LUCOMP,KC2
	Real AMSP,AMK2,Z,FC,A
      EXTERNAL draprn,LUCOMP,ULMASS

C...Flavours fixed when calling from PYREMN or LQQBEV
      IF(IFLRO.EQ.0) GOTO 20

C...Split target remnant qqqQ -> qqQ + q or qqqQbar -> qQbar + qq
C...Q (Qbar) is the partner to the struck sea quark
C...qqq are the nucleon valence quarks from which a quark q or a
C...diquark qq is chosen at random to form a jet system with the
C...scattered sea antiquark or quark, respectively, the other parton
C...forms a baryon qqQ or meson qQbar, respectively.
   10 IFLQ=INT(1.+LST(22)/3.+draprn())
      IF(IFLQ.EQ.LST(22)) THEN
         IFLQQ=2101
         IF(draprn().GT.SNGL(PARL(4))) IFLQQ =2103
      ELSE
         IFLQQ=1000*IFLQ+100*IFLQ+3
      ENDIF
      IFLQ=3-IFLQ

C...Choose flavour of hadron and parton for jet system
      IF(IFLRO.GT.0) THEN
         CALL LUKFDI(IFLQQ,IFLRO,IDUM,K2)
         IF(K2.EQ.0) GOTO 10
         IFLR=IFLQ
      ELSE
         CALL LUKFDI(IFLQ,IFLRO,IDUM,K2)
         IF(K2.EQ.0) GOTO 10
         IFLR=IFLQQ
      ENDIF

C...Entry for use from PYSSPA, flavours given, choose E-p fraction
   20 KSP=IFLR
C...Split energy-momentum of target remnant according to functions P(z)
C...z=E-pz fraction for qq (q) forming jet-system with struck Q (Qbar)
C...1-z=E-pz fraction for qQbar (qqQ) hadron
C...mq=mass of (light) parton remnant q (qq) in jet system
C...mQ=mass of produced (heavy flavour) hadron
      AMSP=ULMASS(KSP)
      AMK2=ULMASS(K2)
C...Old Lepto treatment
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=1/3 since a=1 fixed
      Z=1.-SQRT(draprn())
C...Flip if baryon produced
      KC2=IABS(LUCOMP(K2))
      IF(KC2.GE.301.AND.KC2.LE.400) Z=1.-Z
      IF(LST(14).EQ.2) THEN
C...Update of LEPTO treatment
C...P(z)=(a+1)(1-z)**a with <z>=1/(a+2)=mq/(mq+mQ) --> a=a(mq,mQ)
         A=(AMSP+AMK2)/AMSP - 2.
         Z=draprn()**(1./(A+1.))
      ELSEIF(LST(14).EQ.3) THEN
C...Using Peterson fragmentation function
C...P(z)=N/(z(1-1/z-c/(1-z))**2)  where c=(mq/mQ)**2  (FC=-c)
         FC=-(AMSP/AMK2)**2
   30    Z=draprn()
         IF(-4.*FC*Z*(1.-Z)**2.LT.draprn()*((1.-Z)**2-FC*Z)**2) GOTO 30
      ENDIF
      LST(27)=1
      K2A=IABS(K2)
      IF((K2A.GE.1.AND.K2A.LE.8).OR.K2A.EQ.21.OR.LUCOMP(K2A).EQ.90)
     +LST(27)=2
c      write(6,*) ' LREMH : ',Z,IFLR,K2,AMSP,AMK2
      RETURN
      END
