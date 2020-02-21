C **********************************************************************
	PROGRAM three_Var_th
C
C  COMPUTES THE LEADING CONTRIBUTION IN THE PION CLOUD MODEL TO
C  PROTON PRODUCTION FROM SU(2) MESON-BARYON CONFIGURATIONS OF
C  THE NUCLEON; AS SUCH, THE MAIN MECHANISM IS EXCHANGE OF THE
C  PION AND RHO
C
C  MAINLY, THESE ARE FOR THE REACTION e + n --> p + e' + X, AS
C  MIGHT BE MEASURED AT ``BoNuS''
C
C  MODIFIED FROM PREVIOUS VERSIONS TO INTEGRATE 2 D.O.F. IN THE
C  3D PARAMETER SPACE {x, |k|, theta_h}
C
C  WRITTEN: T. Hobbs (MAY 18, 2014)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nth,ith,ik,FLAG,typ
	PARAMETER (nth=100)
	EXTERNAL fypiN, f_rhoN
	REAL*8  fypiN, f_rhoN
	REAL*8  x,xmin,xmax,xint,k,kmax,kmin,kint,L
	REAL*8  the_h,thmin,thmax,thint,y,kTk
	REAL*8  pi,theta_e,E,mN,Q2,phi_k,rho_PS
	REAL*8  thpt(nth),F2pi0,pi_int
	REAL*8  F2pi_k,F2pi_GRVT,xVpiT,xSpiT
	REAL*8  Ld,F2piK(nth) 
C***********************************************************************
!---------------------------------------------------------------------
!WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
      L = 1.18D0      ! COV DIPOLE: HSS NORM
!      L = 1.71D0      ! IMF DIPOLE: HSS NORM

!      L = 1.63D0       ! HSS +
!      L = 1.56D0      ! HSS CENT. VALUE
!      L = 1.48D0      ! HSS -

!      L = 1.4D0      ! WM 1993 PRD
!...AND USE AN INDEPENDENT PARAMETER FOR THE DELTA:
      Ld = 0.63D0     ! COV DIPOLE: HSS NORM
!      Ld = 1.24D0     ! IMF DIPOLE: HSS NORM

!      Ld = 1.46D0     ! HSS +
!      Ld = 1.39D0     ! HSS CENT. VALUE
!      Ld = 1.32D0     ! HSS -
!_____________________________________________________________________________
!  KINNI
C***********************************************************************
!HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
!      typ = 1 ! DIPOLE FORM FACTOR
!      typ = 2 ! EXPONENTIAL FORM FACTOR
      typ = 3 ! COV. DIPOLE FORM FACTOR
C***********************************************************************
!THIS CHOOSES AMONG THE VARIOUS POSSIBLE COMBINATIONS OF SPLITTING FUNCTION/PDF
      FLAG = 0
!THE FLAGS TOGGLE AMONG DISSOCIATION MODES AS: 
!     FLAG = 0  --- THE PION CONTRIBUTION  | J = 0 + 1/2
!     FLAG = 1  --- THE RHO CONTRIBUTION   | J = 1 + 1/2
!---------------------------------------------------------------------
C***********************************************************************
!EXPERIMENTAL INPUTS AS GIVEN BY THIA KEPPEL (11.12.2013) ------------
       pi = 4.D0*DATAN(1.D0)
       E = 11.D0               !ELECTRON BEAM ENERGY, [GeV]
       theta_e = 35.D0         !ELECTRON SCATTERING ANGLE, DEGREES
!---------------------------------------------------------------------
C***********************************************************************
! WE NOW COMPUTE CONVOLUTIONS AND DISTRIBUTIONS OVER THE RELEVANT theta_h
	ith = 1
	thmin = 0.D0 /180.D0 * pi       !INTEGRATE FULLY OVER pi
	thmax = 180.D0 /180.D0 * pi     !BOUNDS IN HAD. PROD. ANGLE
	thint = (thmax-thmin)/nth
!***************************************
	F2pi0 = 0.D0
	DO   the_h=thmin,thmax+thint/10,thint

	     thpt(ith) = the_h
!---------------------------------------------------------------------
!  CALCULATE USEFUL KINEMATICAL VARIABLES:
       mN = 0.93891897D0          !TARGET NEUTRON MASS, GeV
!---------------------------------------------------------------------
	     F2piK(ith) = 0.D0
! DEFINE THE BOUNDS OF x INTEGRATION
	  ix = 0
	  xmin = 0.05D0
	     xmax = 0.6D0                 !BOUNDS IN ACCESSIBLE x
	     xint = (xmax-xmin)/100.D0    !(FROM THIA!)
	  DO x=xmin,xmax,xint


! DEFINE THE BOUNDS OF THE |k| INTEGRAL -- k1 TO k2
	  ik = 0
	  kmin = 0.06D0
	     kmax = 0.25D0    !SHARP k BOUNDS (FROM THIA!)
	     kint = (kmax-kmin)/100.D0
	  DO k=kmin,kmax,kint

       Q2 = 2.D0*mN*x*E           ![GeV]**2
     & * (1.D0 - 1.D0/( (2.D0*E/(x*mN)) 
     & * DSIN((theta_e/180.D0)*pi)**2.D0 + 1.D0 ) )

       y = (k/mN) * DCOS(the_h)
     &   + (1.D0/mN)*(mN - DSQRT(mN**2.D0 + k**2.D0))
       kTk = k * DSIN(the_h)
       phi_k = DATAN(k/mN)

       rho_PS = (2.D0/mN) * k**2.D0 * DSIN(the_h)  !PHASE SPACE FACTOR
     &        * (1.D0 - DSIN(phi_k)*DCOS(the_h))   !(JACOBIAN)

       CALL GRV (x/y,Q2,xVpiT,xSpiT)
         F2pi_GRVT = (5.D0/9.D0) * (xVpiT + 2.D0 * xSpiT)

!  ------------ INTEGRATED, OVER FINITE y ------------------------
          IF (FLAG.EQ.0) THEN
            pi_int = rho_PS*fypiN(y,kTk,L,typ,0) * F2pi_GRVT
                                                !PION CONTRIBUTION
          ELSE IF (FLAG.EQ.1) THEN
            pi_int = rho_PS*f_rhoN(y,kTk,L,typ,0) * F2pi_GRVT
                                                !RHO CONTRIBUTION
          ENDIF
C________________________________________________________________________________________
C________________________________________________________________________________________
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FOR THE |k| INTEGRAL
	     IF (ik.EQ.0) THEN
	      F2pi_k = F2pi_k + pi_int
	    ELSE IF (ik/2*2.NE.ik) THEN
	      F2pi_k = F2pi_k + 4.D0*pi_int
	    ELSE IF (ik/2*2.EQ.ik) THEN
	      F2pi_k = F2pi_k + 2.D0*pi_int
	     ENDIF

 	    ik = ik + 1
	  ENDDO
	    F2pi_k = (kint/3.D0) * F2pi_k


!----- AND AGAIN FOR THE x INTEGRAL
	     IF (ix.EQ.0) THEN
	      F2piK(ith) = F2piK(ith) + F2pi_k
	    ELSE IF (ix/2*2.NE.ix) THEN
	      F2piK(ith) = F2piK(ith) + 4.D0*F2pi_k
	    ELSE IF (ix/2*2.EQ.ix) THEN
	      F2piK(ith) = F2piK(ith) + 2.D0*F2pi_k
	     ENDIF

 	    ix = ix + 1
	  ENDDO
	    F2piK(ith) = (xint/3.D0) * F2piK(ith)

        print*, the_h, y, F2piK(ith)
!_____________________________________________________________________________
!_____________________________________________________________________________

	  IF (ith/2*2.NE.ith) THEN
	     F2pi0 = F2pi0 + 4.D0*F2piK(ith)
	  ELSE IF (ith/2*2.EQ.ith) THEN
	     F2pi0 = F2pi0 + 2.D0*F2piK(ith)
	  ENDIF
	  ith = ith + 1
	ENDDO
	F2pi0 = (thint/3.D0) * F2pi0

	print*, 'First moment of F2piK(the_h)=',F2pi0
!________________________________________________________________________________________
!________________________________________________________________________________________
C...WRITE DATA TO FILE
!DASSI
          IF (FLAG.EQ.0) THEN
	OPEN (11,FILE='F2pith_INT-x-k.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          ELSE IF (FLAG.EQ.1) THEN
	OPEN (11,FILE='F2rhoth_INT-x-k.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          ENDIF
	  DO ith=2,nth
	  WRITE (11,*) thpt(ith), F2piK(ith)             !CHARGE EXCHANGE
	ENDDO
	CLOSE (11)

	  PRINT*, 'THE NEW DATA HAVE BEEN WRITTEN!'
	END
