*CMZ :  2.05/00 05/03/97  08.08.53  by  Hannes Jung
*CMZ :  2.04/03 16/01/97  11.07.35  by  Hannes Jung
*CMZ :  2.04/00 10/12/96  11.15.26  by  Hannes Jung
*CMZ :  2.03/03 22/08/96  19.17.42  by  Hannes Jung
*CMZ :  2.01/09 10/02/96  19.20.34  by  Hannes Jung
*CMZ :  2.00/01 20/04/95  17.18.33  by  Hannes Jung
*-- Author :
C*********************************************************************

      SUBROUTINE PYSTFE(KF,X,Q2,XPQ)

C...This is a dummy routine, where the user can introduce an interface
C...to his own external structure function parametrization.
C...Arguments in:
C...KF : 11 for e-, 22 for gamma, 211 for pi+, 2212 for p.
C...    Isospin conjugation for n and charge conjugation for
C...    e+, pi-, pbar and nbar is performed by PYSTFU.
C...X : x value.
C...Q2 : Q^2 value.
C...Arguments out:
C...XPQ(-6:6) : x * f(x,Q^2), with index according to KF code,
C...    except that gluon is placed in 0. Thus XPQ(0) = xg,
C...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar,
C...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar,
C...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar,
C...    XPQ(11) = xe-, XPQ(-11) = xe+, XPQ(22) = xgamma.
C...
C...One such interface, to the Diemos, Ferroni, Longo, Martinelli
C...proton structure functions, already comes with the package. What
C...the user needs here is external files with the three routines
C...FXG160, FXG260 and FXG360 of the authors above, plus the
C...interpolation routine FINT, which is part of the CERN library
C...KERNLIB package. To avoid problems with unresolved external
C...references, the external calls are commented in the current
C...version. To enable this option, remove the C* at the beginning
C...of the relevant lines.
C...
C...Alternatively, the routine can be used as an interface to the
C...structure function evolution program of Tung. This can be achieved
C...by removing C* at the beginning of some of the lines below.
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEEP,RGLUDAT2.
      REAL PMAS,PARF,VCKM
      INTEGER KCHG
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
C      SAVE

*KEEP,RGPYPARS.
      INTEGER IPYPAR
      PARAMETER (IPYPAR=200)
      REAL PARP
      INTEGER MSTP
      COMMON/PYPARS/MSTP(IPYPAR),PARP(IPYPAR)
C      SAVE

*KEND.

      DIMENSION XPQ(-6:6)
      CHARACTER CHDFLM(9)*5,HEADER*40
      DATA CHDFLM/'UPVAL','DOVAL','GLUON','QBAR ','UBAR ','SBAR ',
     +'CBAR ','BBAR ','TBAR '/
      DATA HEADER/'Tung evolution package has been invoked'/
      DATA INIT/0/
      IF(MSTP(51).EQ.0.AND.KF.EQ.2212) THEN
         IF(INIT.EQ.0) WRITE(6,*) ' SIMPLE SCALING FUNCTION USED'
         INIT=1
C...Simple scaling structure functions, valence quarks and gluon only.
         XPQ(0)=3.*(1.-X)**5
         XPQ(3)=0.125*(1.-X)**7
         XPQ(2)=70./32 *(1.-X)**3 *SQRT(X)
CCCC    XPQ(2)=(1.-X)**3*(1.274+.589*(1.-X)-1.675*(1.-X)**2)
         XPQ(1)=2.*XPQ(2)
         XPQ(-1)=XPQ(3)
         XPQ(-2)=XPQ(3)
         XPQ(-3)=XPQ(3)
c        RETURN
c      ENDIF

C...Proton structure function evolution from Wu-Ki Tung: parton
C...distribution functions incorporating heavy quark mass effects.
C...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.
      ELSE
         IF(INIT.EQ.0) THEN
            I1=0
            IF(MSTP(52).EQ.4) I1=1
            IHDRN=1
            NU=MSTP(53)
            I2=MSTP(51)
            IF(MSTP(51).GE.11) I2=MSTP(51)-3
            I3=0
            IF(MSTP(52).EQ.3) I3=1

C...Convert to Lambda in CWZ scheme (approximately linear relation).
            ALAM=0.75*PARP(1)
            TPMS=PMAS(6,1)
            QINI=PARP(52)
            QMAX=PARP(53)
            XMIN=PARP(54)

C...Initialize evolution (perform calculation or read results from
C...file).
C...Remove C* on following two lines to enable Tung initialization.
C*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,
C*   &    I2,I3,IRET,IRR)
            INIT=1
         ENDIF

C...Put into output array.
         Q=SQRT(Q2)
         DO 10 I=-6,6
            FIXQ=0.
C...Remove C* on following line to enable structure function call.
C*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR))
   10    XPQ(I)=X*FIXQ

C...Change order of u and d quarks from Tung to PYTHIA convention.
         XPS=XPQ(1)
         XPQ(1)=XPQ(2)
         XPQ(2)=XPS
         XPS=XPQ(-1)
         XPQ(-1)=XPQ(-2)
         XPQ(-2)=XPS
      ENDIF

      RETURN
      END
