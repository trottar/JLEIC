C **********************************************************************
	PROGRAM three_Var_theta
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
C  MODIFIED, MAY 8, 2015
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nth,ith,ik,FLAG,typ,dis
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
!      L = 1.18D0      !COV DIPOLE: HSS NORM
!      L = 1.71D0      !IMF DIPOLE: HSS NORM

!      L = 1.63D0       !HSS +
      L = 1.56D0      !HSS CENT. VALUE
!      L = 1.48D0      !HSS -

!      L = 1.4D0      !WM 1993 PRD
!...AND USE AN INDEPENDENT PARAMETER FOR THE DELTA:
!      Ld = 0.63D0     !COV DIPOLE: HSS NORM
!      Ld = 1.24D0     !IMF DIPOLE: HSS NORM

!      Ld = 1.46D0     !HSS +
      Ld = 1.39D0     !HSS CENT. VALUE
!      Ld = 1.32D0     !HSS -
!_____________________________________________________________________________
!  KINNI
C***********************************************************************
!HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
!      typ = 1 !DIPOLE FORM FACTOR
      typ = 2 !EXPONENTIAL FORM FACTOR
!      typ = 3 !COV. DIPOLE FORM FACTOR
C***********************************************************************
!THE PARAMETER 'dis' SELECTS THE CHARGE CONFIG. OF THE INTERM. STATE
!       dis = 0 !CHARGE EXCHANGE
       dis = 1 !NEUTRAL EXCHANGE
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
            pi_int = rho_PS*fypiN(y,kTk,L,typ,dis) * F2pi_GRVT
                                                !PION CONTRIBUTION
          ELSE IF (FLAG.EQ.1) THEN
            pi_int = rho_PS*f_rhoN(y,kTk,L,typ,dis) * F2pi_GRVT
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

        print*, the_h*180/3/14, k, F2piK(ith)
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
	OPEN (11,FILE='BONUS-3Var/F2pith_INT-x-k.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          ELSE IF (FLAG.EQ.1) THEN
	OPEN (11,FILE='BONUS-3Var/F2rhoth_INT-x-k.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          ENDIF
	  DO ith=2,nth
	  WRITE (11,*) thpt(ith)*(180.D0/pi), F2piK(ith)
	ENDDO
	CLOSE (11)

	  PRINT*, 'THE NEW DATA HAVE BEEN WRITTEN!'
	END
C ***********************************************************************
C ***********************************************************************
	FUNCTION fypiN (y,kT,L,typ,dis)
C
C  FUNCTION GIVING f(y) FOR N-PION VERTEX, WHERE
C    y IS THE IMF MOMENTUM OF THE INTERMEDIATE MESON.
C
C  WRITTEN: W. Melnitchouk (1999)
C  MODIFIED: T. HOBBS (2013)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ,dis
	REAL*8  ss,kT,kT2,SpiN
	REAL*8  y,fypiN,t
	REAL*8  pi,mN,mpi,mP,g_piNN,gg,FF,L,sM

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> N + pi^+, or N --> P + pi^- 
	mpi = 0.13957018D0    !THE MASS OF THE pi^+/-
	mP = mN               !THE MASS OF THE NUCLEON
       g_piNN = DSQRT (14.40D0 * 4*pi)    ! g_{pi NN}
	gg = 2.D0 * g_piNN**2 / (16.D0 * pi**2)  !2 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> P + pi^0
	mpi = 0.13957018D0    !THE MASS OF THE pi^0
	mP = mN               !THE MASS OF THE NUCLEON
        g_piNN = DSQRT (14.40D0 * 4*pi)    ! g_{pi NN}
	gg = g_piNN**2 / (16.D0 * pi**2)  !1 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF
	IF (y.LE.0.D0 .OR. y.GE.1.D0) THEN
	  fypiN = 0.D0
	  RETURN
	ENDIF

           kT2 = kT**2
           SpiN = (kT2 + mpi**2)/y + (kT2 + mP**2)/(1.D0-y)

          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + SpiN))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + SpiN))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SpiN)/(L**2) )          ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mpi**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN                       ! DIPOLE -- s-channel Lambda exchange
            sM = (kT2 + (1.D0+y)*mpi**2.D0)/y + (kT2 
     &                         + y*mP**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mP**4)/(L**4 + sM**2) 
          ENDIF

        ss = ( kT2 + (mP - (1.D0-y) * mN)**2.D0) / (1.D0-y)
     &       / ( (1.D0-y)*(SpiN - mN**2.D0) )**2.D0 * FF**2.D0

!        ss = ( kT2 + (mP - (1.D0-y) * mN)**2.D0) / (1.D0-y)  !NO FORM FACTOR!
!     &       / ( (1.D0-y)*(SpiN - mN**2.D0) )**2.D0

        fypiN =  gg * (1.D0-y) / y * ss
	RETURN   
	END
C ***************************************************************************
        FUNCTION f_rhoN (y,kT,L,typ,dis)
C  Function giving numerical value of f(y) for N-rho-N
C  OUTPUT IS THE NEUTRON ---> rho^- + PROTON SPLITTING FUNCTION
C  By: T. Hobbs on NOV 12, 2013
C  taken from notes "spin-1, m_B /= m_N," May 17, 2012
C ***************************************************************************
        IMPLICIT NONE
        INTEGER typ,dis
        REAL*8  kT,kT2,SRoN,P_k,pl_k,P_p,t,sv,st,si,ss
        REAL*8  y,f_rhoN,sM
        REAL*8  pi,mN,mrho,mP,g_RoNN,f_RoNN,gg,fff,fg,FF,L
     
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!**** THE DISSOCIATION N --> rho^- + PROTON
	mrho = 0.7754D0  !THE MASS OF THE rho^-
	mP = mN  !THE MASS OF THE PROTON
!*** WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
	g_RoNN = DSQRT (2.D0 * 0.55D0 * 4.D0 * pi)  !rho-N-N (Hohler and Pieteranin)
	f_RoNN = 6.1D0 * g_RoNN
	gg = g_RoNN**2 / (16.D0 * pi**2)
	fff = f_RoNN**2 / (16.D0 * pi**2)
	fg = f_RoNN*g_RoNN / (16.D0 * pi**2) !FACTOR OF 2: ISOSPIN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSEIF (dis.EQ.1) THEN
!**** THE DISSOCIATION P --> rho^0 + P
	mrho = 0.7754D0  !THE MASS OF THE rho^0
	mP = mN  !THE MASS OF THE PROTON
!*** WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
	g_RoNN = DSQRT (0.55D0 * 4.D0 * pi)  !rho-N-N (Hohler and Pieteranin)
	f_RoNN = 6.1D0 * g_RoNN
	gg = g_RoNN**2 / (16.D0 * pi**2)
	fff = f_RoNN**2 / (16.D0 * pi**2)
	fg = f_RoNN*g_RoNN / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF

        IF (y.LE.0.D0 .OR. y.GE.1.D0) THEN
          f_rhoN = 0.D0
          RETURN
        ENDIF

           kT2 = kT**2
           SRoN = (kT2 + mrho**2)/y + (kT2 + mP**2)/(1.D0-y)
     
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))       ! MONOPOLE
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))**2    ! DIPOLE
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SRoN)/(L**2) )          ! EXPONENTIAL
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mrho**2) / (L**2 - t))**2     ! COV. DIPOLE
          ELSE IF (typ.EQ.4) THEN                       ! DIPOLE -- s-CHANNEL LAMBDA EXCHANGE
            sM = (kT2 + (1.D0+y)*mrho**2.D0)/y + (kT2 
     &                         + y*mP**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mP**4)/(L**4 + sM**2) 
          ENDIF

C...TOPT WITH P[alpha] - p[alpha] DERIVATIVE COUPLING
          P_k = (mrho**2 + y**2*mN**2 + kT2)/2.D0/y
          P_p = (mP**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
          pl_k = (mP**2+kT2)*y/2.D0/(1.D0-y) 
     &          + (mrho**2+kT2)*(1.D0-y)/2.D0/y + kT2  

          sv = -6.D0*mN*mP + 4.D0*P_k*pl_k/mrho**2 + 2.D0*P_p

          st = -(P_p)**2 + P_p*(mP+mN)**2 - mP*mN*(mP**2+mN**2+mP*mN)
     &       + 1.D0/(2.D0*mrho**2) * ( (P_p - mP*mN)*(P_k-pl_k)**2
     &       - 2.D0*(P_k-pl_k)*(mP**2*P_k - mN**2*pl_k)
     &       + 2.D0*P_k*pl_k*(2.D0*P_p-mP**2-mN**2) ) 

          si = -4.D0*(mP+mN)*(mP*mN - P_p) 
     &         -2.D0*(mP*P_k**2 - (mP+mN)*P_k*pl_k + mN*pl_k**2)/mrho**2
          
          IF (kt.gt.0.0 .and. (sv.lt.0.0 .or. st.lt.0.0)) THEN
            PRINT *,'CS1 -- ##### kT,y,sv,st =',kt,y,sv,st
            STOP 
          ENDIF
        
          ss = (gg*sv + fff*st/mN**2 + fg*si/mN)    !THE FULL
     &       / ( y*(SRoN - mN**2) )**2 * FF**2      !EXPRESSION (UN-INT.)**
!__________________________________________________________________________
!************ TO STUDY SEPARATE TERMS' CONTRIBUTIONS!!!! ******************
!          ss = gg*sv                                         
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !VECTOR
!          ss = fg*si/mN                                      
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !INTERFERENCE
!          ss = fff*st/mN**2                                  
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !TENSOR
!__________________________________________________________________________
        f_rhoN = y / (1.D0-y) * ss
        RETURN
        END
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C ***************************************************************************
C ***************************************************************************
        FUNCTION f_RhoDel (y,kT,L,typ,dis)
C  Function giving numerical value of f(y) for N-Rho-Delta
C  OUTPUT IS THE NEUTRON ---> RHO + DELTA SPLITTING FUNCTION
C  By: T. Hobbs on JAN 7, 2014
C ***************************************************************************
        IMPLICIT NONE
        INTEGER typ,dis
        REAL*8  kT,kT2,SRoD,P_k,pl_k,P_p,t
        REAL*8  y,f_RhoDel,sr,ss
        REAL*8  pi,mN,mD,mRo,g_NDS,gg,FF,L,sM
     
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION N --> (RHO^0 + DELTA^0) + (RHO^- + DELTA^+)
        mRo = 0.7754D0   !THE MASS OF THE RHO MESON
        mD = 1.232D0 !Mass OF THE DELTA ISOBAR
!!*** WE USE THE COUPLINGS INFERRED FROM HOLZENKAMP ET AL. ***
        g_NDS = DSQRT (20.448D0 * 4.D0*pi)    ! as N-D*-Sigma*_c
        gg = (2.D0/3.D0) * g_NDS**2 / (16.D0 * pi**2 * mRo**2)
!WE INCLUDE AN OVERALL FACTOR OF 2/3 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** A BLANK SPACE FOR OTHER SEPARATE MODES
        mRo = 0.7754D0   !THE MASS OF THE RHO MESON
        mD = 1.232D0 !Mass OF THE DELTA BARYON
!!*** TERMINATE THE OUTPUT ***
        g_NDS = 0.D0    ! ZERO OUTPUT
        gg = g_NDS**2 / (16.D0 * pi**2 * mRO**2)
!WE INCLUDE AN OVERALL FACTOR OF 1 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF
        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          f_RhoDel = 0.D0
          RETURN
        ENDIF
     
          kT2 = kT**2
          SRoD = (kT2 + mRo**2)/y + (kT2 + mD**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoD))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SRoD)/(L**2) )        ! expon
          ELSE IF (typ.EQ.3) THEN
            t = ( -kT2 - (1.D0-y)*(mRo**2 - y*mN**2) ) / y
            FF = ((L**2 - mD**2) / (L**2 - t))**2      ! t-channel dip
          ELSE IF (typ.EQ.4) THEN
	  sM = (kT2 + (2.D0-y)*mRo**2.D0)/(1.D0-y) + (kT2 + (1.D0-y)*mD**2.D0)
     &       /y + mN**2.D0
          FF = (L**4 + mD**4)/(L**4 + sM**2)  !DIPOLE -- s-channel Lambda exchange
          ENDIF

C...TOPT with P[alpha] - p[alpha] derivative coupling
          P_k = (mRo**2 + y**2*mN**2 + kT2)/2.D0/y
          P_p = (mD**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
          pl_k = (mD**2+kT2)*y/2.D0/(1.D0-y) 
     &          + (mRo**2+kT2)*(1.D0-y)/2.D0/y + kT2  

         sr = -4.D0*mN*mD/3.D0*(2.D0*mD**2.D0+mN*mD+2.D0*mN**2.D0)
     &   -4.D0*mN*mD/(3.D0*mRo**2.D0)*(P_k-pl_k)**2.D0
     &   -4.D0/(3.D0*mRo**2.D0)*(mD**2.D0*P_k**2.D0+mN**2.D0*pl_k**2.D0)
     &   +4.D0*P_p/3.D0*(2.D0*mD**2.D0+4.D0*mN*mD+mN**2.D0)
     &   +4.D0*P_p/(3.D0*mRo**2.D0)*pl_k**2.D0*(1.D0-mN**2.D0/mD**2.D0)
     &   -4.D0*P_p**2.D0*(1.D0-2.D0*P_k*pl_k/(3.D0*mRo**2.D0*mD**2.D0)
     &       -P_p/(3.D0*mD**2.D0))
          
          if (kt.gt.0.0 .and. (sr.lt.0.0)) then
            print *,'CS2 -- ##### kT,y,sr =',kt,y,sr
            stop
          endif
        
          ss = sr
     &       /((1.D0-y)*(SRoD - mN**2))**2 * FF**2

        f_RhoDel = gg * (1.D0-y) / y * ss

        RETURN
        END
C ***************************************************************************
C ***********************************************************************
	FUNCTION fypiD (y,kT,L,typ,dis)
C
C  Function giving numerical value of f(y) for the SU(4) analogue of the pi-Delta
C  interaction, as usual y is IMF momentum fraction of the charmed BARYON.
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ,dis
	REAL*8  ss,kT,kT2,SpiD
	REAL*8  y,yR,fypiD,t
	REAL*8  pi,mN,mD,mpi,g_DLcN,gg,FF,L,sM

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION N ---> (pi^0 + DELTA^0) + (pi^- + DELTA^-)
	mD = 1.232D0  !THE MASS OF THE DELTA BARYON
	mpi = 0.13957018D0 !Mass OF THE PION 
!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
	g_DLcN = DSQRT (0.2237D0 * 4*pi)    ! as N-pi-Del from Holzenkamp et al.
	gg =  (2.D0/3.D0) * g_DLcN**2 / (16.D0 * pi**2)   !2/3 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE NULL DISSOCIATION
	mD = 1.232D0  !THE MASS OF THE DELTA
	mpi = 0.13957018D0 !Mass OF THE PION
!!*** WE USE ZERO COUPLINGS TO RETURN A NULL OUTPUT
	g_DLcN = 0.D0    
	gg = g_DLcN**2 / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF

	IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
	  fypiD = 0.D0
	  RETURN
	ENDIF

	  kT2 = kT**2
	  SpiD = (kT2 + mpi**2)/y + (kT2 + mD**2)/(1.D0-y)
!********************************************************************************
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + SpiD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + SpiD))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SpiD)/(L**2) )          ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mpi**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN                       ! DIPOLE -- s-channel Lambda exchange
            sM = (kT2 + (1.D0+y)*mpi**2.D0)/y + (kT2 
     &                         + y*mD**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mD**4)/(L**4 + sM**2) 
          ENDIF

         yR = 1.D0 - y
      ss = ( kT2 + (mD-yR*mN)**2) * ( kT2 + (mD+yR*mN)**2 )**2
     &       / ( 6.D0*mD**2*yR**3) / ( (1.D0-yR)*(SpiD - mN**2) )**2 
     &       * FF**2

	fypiD =  gg/(mpi**2) * (1.D0-yR) / yR * ss

	RETURN   
	END
C ***************************************************************************
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C        THE PION STRUCTURE FUNCTION IS TAKEN FROM THIS SUBROUTINE
C ***************************************************************************
        SUBROUTINE GRV (x,Q2,xVpi,xSpi)
C
C  Subroutine giving x-dependence of parametrization of LO pion valence
C  and sea distribution functions xVpi, xSpi, valid for 0.25 < Q2 < 10^8 GeV^2,
C  and 10^-5 < x < 1.
C
C  Gluck, Reya, Vogt: Z.Phys. C53 (1992) 651 (appendix 1).
C ******************************************************************************
!        IMPLICIT UNDEFINED (A-Z)
!        IMPLICIT NONE (A-Z)
        REAL*8  x,xVpi,a,AA,D,Nv,
     &          xSpi,alpha,as,AAs,Bs,Ds,E,Epr,beta
        REAL*8  Q2,Q02,L,s

        xVpi = 0.D0
        xSpi = 0.D0

	IF (x.LT.1.D-5) x=1.01D-5
	IF (Q2.LT.0.25D0) Q2=0.2501D0

        Q02 = 0.25D0
        L = 0.232D0
        IF (Q2.LE.Q02) RETURN
        s = DLOG( DLOG(Q2/L**2) / DLOG(Q02/L**2) )

C...Valence distribution
        Nv = 0.519D0 + 0.180D0*s - 0.011D0*s**2
        a = 0.499D0 - 0.027D0*s
        AA = 0.381D0 - 0.419D0*s
        D = 0.367D0 + 0.563D0*s
        xVpi = 0.D0
        IF (x.LT.1.D0)
     &	  xVpi = Nv * x**a * (1.D0+AA*DSQRT(x)) * (1.D0-x)**D

C...Sea distribution (SU(3) symmetric)
        alpha = 0.55D0
        as = 2.538D0 - 0.763D0*s
        AAs = -0.748D0
        Bs = 0.313D0 + 0.935D0*s
        Ds = 3.359D0
        E = 4.433D0 + 1.301D0*s
        Epr = 9.30D0 - 0.887D0*s
        beta = 0.56D0
        xSpi = 0.D0
        IF (x.LT.1.D0)
     &	  xSpi = s**alpha / (DLOG(1.D0/x))**as
     &         * (1.D0 + AAs*DSQRT(x) + Bs*x) * (1.D0 - x)**Ds
     &         * DEXP(-E + DSQRT(Epr*s**beta*DLOG(1.D0/x)))

        RETURN
        END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
