C **********************************************************************
	PROGRAM CS_gfor
C
C-----  OUTPUTS THE NUMERICAL VALUES OF MBM HIGH X_F INCLUSIVE CROSS
C-----  SECTIONS FOR THE PURPOSE OF EXP. DATA COMPARISONS
C
C  Written: T. Hobbs (Sept 19, 2013)
C  MADE COMPATIBLE WITH gfortran: TH (June 18, 2014)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,FLAG,typ
	PARAMETER (nx=100)
	EXTERNAL fyLcD, f_DstL, f_DstSig, fypiD
	REAL*8  fypiD, fyLcD, f_DstL, f_DstSig
	REAL*8  x,xmin,xmax,xint,xpt(nx)
	REAL*8  L,pi,s_TOT,sig(nx)
C***********************************************************************
!---------------------------------------------------------------------
!WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
      L = 1.003D0  !FOR LAMBDA    (L = [1.003 +/- 0.008] GeV)
!          L = 1.011 !UPPER BOUND --- STAT.  | -- LAMBDA PRODUCTION
!          L = 0.995 !LOWER BOUND --- STAT.  |
!      L = 1.17D0  !FOR SIGMA+    (L = [1.17 +/- 0.01] GeV)
!      L = 1.24D0  !FOR SIGMA*+   (L = [1.24 +/- 0.02] GeV)
C***********************************************************************
!HERE WE PLACE A GLOBAL FLAG FOR THE SUPPRESSION FACTOR
      typ = 1      !EXPONENTIAL FORM FACTOR: HOLTMANN ET AL. (1996)
!      typ = 2       !EXPONENTIAL FORM FACTOR
C***********************************************************************
!THIS CHOOSES AMONG THE VARIOUS FINAL STATE BARYONS
      FLAG = 0
!THE FLAGS TOGGLE AMONG DIFFERENT PRODUCTION MECHANISMS 
!     FLAG = 0  --- (pp --> Lambda + X)
!     FLAG = 1  --- (pp --> Sigma+ + X)
!     FLAG = 2  --- (pp --> Sigma*+ + X)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! WE NOW COMPUTE CROSS SECTIONS OVER A FULL RANGE OF x
	xmax = 1.D0
	xint = xmax/DBLE(nx)

        DO ix = 0, nx
            x = xint * DBLE(ix)
            xpt(ix) = x

      pi = 4*DATAN(1.D0)
      s_TOT = 19.9D0

          IF (FLAG.EQ.0) THEN
      sig(ix) = s_TOT * (x/pi) * (f_DstL(1.D0-x,L,typ,0) 
     &       + f_DstL(1.D0-x,L,typ,1) 
     &       + fyLcD(x,L,typ,0) + fyLcD(x,L,typ,1)
     &       + f_DstSig(1.D0-x,L,typ,0)
     &       + fypiD(x,L,typ,0) )
!---------------------------------------------------------------
          ELSE IF (FLAG.EQ.1) THEN
      sig(ix) = s_TOT * (x/pi) 
     &       * ( fyLcD(x,L,typ,2) + f_DstL(1.D0-x,L,typ,2) )
!---------------------------------------------------------------
          ELSE IF (FLAG.EQ.2) THEN
      sig(ix) = s_TOT * (x/pi)
     &       * ( fypiD(x,L,typ,1) + f_DstSig(1.D0-x,L,typ,1) )
!---------------------------------------------------------------
          ENDIF
        ENDDO
C*******************************************************************
C...Write to file
          IF (FLAG.EQ.0) THEN
       OPEN (10,FILE='Lambda_MBM.dat',STATUS='UNKNOWN',
     &                                           FORM='FORMATTED')
          ELSE IF (FLAG.EQ.1) THEN
       OPEN (10,FILE='Sigma_MBM.dat',STATUS='UNKNOWN',
     &                                           FORM='FORMATTED')
          ELSE IF (FLAG.EQ.2) THEN
       OPEN (10,FILE='Sigma-st_MBM.dat',STATUS='UNKNOWN',
     &                                           FORM='FORMATTED')
          ENDIF
	  DO ix=1,nx
	  WRITE (10,*) xpt(ix), sig(ix)
	ENDDO
	CLOSE (10)
	     print*, 'THE DATA ARE WRITTEN'
	END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C ***********************************************************************
	FUNCTION fyLcD (y,L,typ,dis)
C
C  Function giving numerical value of f(y) for N-Lambda-K vertex
C    y is l.c. momentum fraction on baryon.
C
C  Modified: T. Hobbs (2012)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER ikT,typ,dis
	REAL*8  ss,ss0,kT,kT2,kTmax,kTint,SLcD
	REAL*8  y,fyLcD,t
	REAL*8  pi,mN,mD0,mLc,g_DLcN,gg,FF,L,sM

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> K^+ + LAMBDA^0 
	mD0 = 0.4937D0  !THE MASS OF THE K+
	mLc = 1.1157D0  !THE MASS OF THE LAMBDA^0
!!*** WE USE THE COUPLINGS INFERRED FROM MUELLER-GROELING ET AL. ***
	g_DLcN = DSQRT (15.56D0 * 4*pi)    ! Mueller-Groeling - PS
!       g_DLcN = DSQRT (14.40D0 * 4*pi)    ! g_{pi NN}
!	g_DLcN = DSQRT (9.D0 * 4*pi)       ! Navarra I
!	g_DLcN = DSQRT (3.61D0 * 4*pi)     ! Navarra II
	gg = g_DLcN**2 / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> K^+ + SIGMA^0
	mD0 = 0.4937D0   !THE MASS OF THE K+
	mLc = 1.1926D0  !THE MASS OF THE SIGMA^0
!!*** WE USE THE COUPLINGS INFERRED FROM M-G ET AL. ***
	g_DLcN = DSQRT (0.576D0 * 4*pi)    ! as N-K-Sigma
	gg = g_DLcN**2 / (16.D0 * pi**2)  !1 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.2) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> K^0 + SIGMA^+
	mD0 = 0.4976D0   !THE MASS OF THE K^0
	mLc = 1.1894D0  !THE MASS OF THE SIGMA^+
!*** WE USE THE COUPLINGS INFERRED FROM M-G ET AL. ***
	g_DLcN = DSQRT (0.576D0 * 4*pi)    ! as N-K-Sigma
	gg = 2.D0*g_DLcN**2 / (16.D0 * pi**2)  !2 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF
	IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
	  fyLcD = 0.D0
	  RETURN
	ENDIF

	ss0 = 0.D0
	kTmax = 10.D0
	kTint = kTmax/1000.D0

        DO ikT = 0, 1000
            kT = kTint * DBLE(ikT)

	  kT2 = kT**2
	  SLcD = (kT2 + mD0**2)/(1.D0-y) + (kT2 + mLc**2)/y

          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + SLcD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = DEXP( (mN**2 - SLcD)/(2*L**2) )        ! expon -- HSS
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SLcD)/(L**2) )        ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mD0**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN                      ! DIPOLE -- s-channel Lambda exchange
            sM = (kT2 + (1.D0+y)*mD0**2.D0)/y + (kT2 
     &                         + y*mLc**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mLc**4)/(L**4 + sM**2) 
          ENDIF

      ss = ( kT2 + (mLc-y*mN)**2) / y
     &       / ( (1.D0-y)*(SLcD - mN**2) )**2 * FF**2 * (2*kT)

	  IF (ikT/2*2.NE.ikT) THEN
	    ss0 = ss0 + 4*ss
	  ELSE IF (ikT/2*2.EQ.ikT) THEN
	    ss0 = ss0 + 2*ss
	  ENDIF

	ENDDO
	fyLcD =  gg * (1.D0-y) / y * (kTint/3) * ss0
	RETURN   
	END
C ***************************************************************************
        FUNCTION f_DstL (y,L,typ,dis)
C  Function giving numerical value of f(y) for N-K*-Lambda
C  OUTPUT IS THE PROTON ---> K* + Lambda SPLITTING FUNCTION
C  By: T. Hobbs on June 6, 2012
C  taken from notes "spin-1, m_B /= m_N," May 17, 2012
C ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SRoN,P_k,pl_k,P_p,t,
     &          sv,st,si,ss,ss0
        REAL*8  y,f_DstL,sM
        REAL*8  pi,mN,mD,mL,g_RoNN,f_RoNN,gg,fff,fg,FF,L
     
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!**** THE DISSOCIATION P --> K*+ + LAMBDA^0
	mD = 0.8917D0  !THE MASS OF THE K*+
	mL = 1.1157D0  !THE MASS OF THE LAMBDA^0
!*** WE USE THE COUPLINGS INFERRED FROM M-G ET AL. ***
	g_RoNN = DSQRT (2.503D0 * 4.D0 * pi)  !N-K*-Lambda
	f_RoNN = 3.259D0 * g_RoNN
	gg = g_RoNN**2 / (16.D0 * pi**2)
	fff = f_RoNN**2 / (16.D0 * pi**2)
	fg = f_RoNN*g_RoNN / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!**** THE DISSOCIATION P --> K*+ + SIGMA^0
	mD = 0.8917D0  !THE MASS OF THE K*+
	mL = 1.1926D0  !THE MASS OF THE SIGMA^0
!*** WE USE THE COUPLINGS INFERRED FROM M-G ET AL. ***
	g_RoNN = DSQRT (0.841D0 * 4.D0 * pi)  !N-K*-Sigma
	f_RoNN = -2.42D0 * g_RoNN !NOTE THE RELATIVE MINUS SIGN!!
	gg = g_RoNN**2 / (16.D0 * pi**2)
	fff = f_RoNN**2 / (16.D0 * pi**2)
	fg = f_RoNN*g_RoNN / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.2) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!**** THE DISSOCIATION P --> K*0 + SIGMA^+ 
	mD = 0.8958D0   !THE MASS OF THE K*0
	mL = 1.1894D0  !THE MASS OF THE SIGMA^+
!*** WE USE THE COUPLINGS FROM MUELLER -- K*-SIG  ---
	g_RoNN = DSQRT (0.841D0 * 4.D0 * pi * 2.D0) !2 ISOSPIN
	f_RoNN = -2.42D0 * g_RoNN
	gg = g_RoNN**2 / (16.D0 * pi**2)
	fff = f_RoNN**2 / (16.D0 * pi**2)
	fg = f_RoNN*g_RoNN / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF

        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          f_DstL = 0.D0
          RETURN
        ENDIF
     
        ss0 = 0.D0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  

        DO ikT = 0, 1000
            kT = kTint * DBLE(ikT)

          kT2 = kT**2
          SRoN = (kT2 + mD**2)/y + (kT2 + mL**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = DEXP( (mN**2 - SRoN)/(2*L**2) )        ! expon -- HSS
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SRoN)/(L**2) )        ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mD**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN                      ! DIPOLE -- s-channel Lambda exchange
            sM = (kT2 + (1.D0+y)*mD**2.D0)/y + (kT2 
     &                         + y*mL**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mL**4)/(L**4 + sM**2) 
          ENDIF

C...TOPT with P[alpha] - p[alpha] derivative coupling
          P_k = (mD**2 + y**2*mN**2 + kT2)/2.D0/y
          P_p = (mL**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
          pl_k = (mL**2+kT2)*y/2.D0/(1.D0-y) 
     &          + (mD**2+kT2)*(1.D0-y)/2.D0/y + kT2  

          sv = -6.D0*mN*mL + 4.D0*P_k*pl_k/mD**2 + 2.D0*P_p

          st = -(P_p)**2 + P_p*(mL+mN)**2 - mL*mN*(mL**2+mN**2+mL*mN)
     &       + 1.D0/(2.D0*mD**2) * ( (P_p - mL*mN)*(P_k-pl_k)**2
     &       - 2.D0*(P_k-pl_k)*(mL**2*P_k - mN**2*pl_k)
     &       + 2.D0*P_k*pl_k*(2.D0*P_p-mL**2-mN**2) ) 

          si = -4.D0*(mL+mN)*(mL*mN - P_p) 
     &         -2.D0*(mL*P_k**2 - (mL+mN)*P_k*pl_k + mN*pl_k**2)/mD**2
          
          if (kt.gt.0.0 .and. (sv.lt.0.0 .or. st.lt.0.0)) then
            print *,'CS1 -- ##### kT,y,sv,st =',kt,y,sv,st
            stop
          endif
        
          ss = (gg*sv + fff*st/mN**2 + fg*si/mN)             !FULL
     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)


          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
        ENDDO
        f_DstL = y / (1.D0-y) * (kTint/3.D0) * ss0 
        RETURN
        END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        FUNCTION f_DstSig (y,L,typ,dis)
C  Function giving numerical value of f(y) for N-D*-Sigmaa_c^+
C  OUTPUT IS THE PROTON ---> K* + Sigma* SPLITTING FUNCTION
C  By: T. Hobbs; May, 2014
C  taken from notes WM, Japanese Proceedings, 1993
C ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SRoN,P_k,pl_k,P_p,t
        REAL*8  y,f_DstSig,sr,ss,ss0
        REAL*8  pi,mN,mD,mS,g_NDS,gg,FF,L,sM
     
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> K*+ + SIGMA* 
        mD = 0.8917D0  !THE MASS OF THE K*+
        mS = 1.3837D0  !MASS OF THE NEUTRAL Sigma*
!!*** WE USE THE COUPLINGS INFERRED FROM HOLZENKAMP ET AL. ***
        g_NDS = DSQRT (3.408D0 * 4*pi)    ! as N-K*-Sigma*
        gg = g_NDS**2 / (16.D0 * pi**2 * mD**2)
!WE INCLUDE AN OVERALL FACTOR OF 1 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> K* + SIGMA*+ 
        mD = 0.8958D0   !THE MASS OF THE \bar{K}*0
        mS = 1.3828D0   !Mass of the Sigma*^+
!!*** WE USE THE COUPLINGS INFERRED FROM HOLZENKAMP ET AL. ***
        g_NDS = DSQRT (3.408D0 * 4*pi * 2.D0)    ! as N-K*-Sigma*
        gg = g_NDS**2 / (16.D0 * pi**2 * mD**2)
!WE INCLUDE AN OVERALL FACTOR OF 2 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF
        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          f_DstSig = 0.D0
          RETURN
        ENDIF
     
        ss0 = 0.D0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  

        DO ikT = 0, 1000
            kT = kTint * DBLE(ikT)
          kT2 = kT**2
          
          SRoN = (kT2 + mD**2)/y + (kT2 + mS**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = DEXP( (mN**2 - SRoN)/(2*L**2) )        ! expon -- HSS
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SRoN)/(L**2) )        ! expon
          ELSE IF (typ.EQ.3) THEN
            t = ( -kT2 - (1.D0-y)*(mD**2 - y*mN**2) ) / y
            FF = ((L**2 - mS**2) / (L**2 - t))**2      ! t-channel dip
          ELSE IF (typ.EQ.4) THEN
      sM = (kT2 + (2.D0-y)*mD**2.D0)/(1.D0-y)
     &   + (kT2 + (1.D0-y)*mS**2.D0)/y + mN**2.D0
          FF = (L**4 + mS**4)/(L**4 + sM**2)  !DIPOLE -- s-channel Lambda exchange
          ENDIF

C...TOPT with P[alpha] - p[alpha] derivative coupling
          P_k = (mD**2 + y**2*mN**2 + kT2)/2.D0/y
          P_p = (mS**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
          pl_k = (mS**2+kT2)*y/2.D0/(1.D0-y) 
     &          + (mD**2+kT2)*(1.D0-y)/2.D0/y + kT2  

         sr = -4.D0*mN*mS/3.D0*(2.D0*mS**2.D0+mN*mS+2.D0*mN**2.D0)
     &    -4.D0*mN*mS/(3.D0*mD**2.D0)*(P_k-pl_k)**2.D0
     &    -4.D0/(3.D0*mD**2.D0)*(mS**2.D0*P_k**2.D0+mN**2.D0*pl_k**2.D0)
     &    +4.D0*P_p/3.D0*(2.D0*mS**2.D0+4.D0*mN*mS+mN**2.D0)
     &    +4.D0*P_p/(3.D0*mD**2.D0)*pl_k**2.D0*(1.D0-mN**2.D0/mS**2.D0)
     &    -4.D0*P_p**2.D0*(1.D0-2.D0*P_k*pl_k/(3.D0*mD**2.D0*mS**2.D0)
     &       -P_p/(3.D0*mS**2.D0))
          
          if (kt.gt.0.0 .and. (sr.lt.0.0)) then
            print *,'CS2 -- ##### kT,y,sr =',kt,y,sr
            stop
          endif
        
          ss = sr
     &        /((1.D0-y)*(SRoN - mN**2))**2 * FF**2  * (2.D0*kT)

          IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF

        ENDDO
        f_DstSig = gg * (1.D0-y) / y * (kTint/3.D0) * ss0 

        RETURN
        END
C ***************************************************************************
C ***********************************************************************
	FUNCTION fypiD (y,L,typ,dis)
C
C  Function giving numerical value of f(y) for the SU(3) analogue of the pi-Delta
C  interaction, as usual y is IMF momentum fraction of the strange BARYON.
C ***********************************************************************
	IMPLICIT NONE
	INTEGER ikT,typ,dis
	REAL*8  ss,ss0,kT,kT2,kTmax,kTint,SLcD
	REAL*8  y,fypiD,t,sM
	REAL*8  pi,mN,mD0,mLc,g_DLcN,gg,FF,L

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0 !masses in GeV!!
          IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> SIGMA* + K+
	mD0 = 0.4937D0   !THE MASS OF THE K+
	mLc = 1.3837D0   !Mass of the NEUTRAL Sigma*
	g_DLcN = DSQRT (0.0372D0 * 4*pi)    ! as N-Y*-K from Holzenkamp et al.
	gg =  g_DLcN**2 / (16.D0 * pi**2)   !1 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!**** THE DISSOCIATION P --> SIGMA*+ + \bar{K}^0 
	mD0 = 0.4976D0   !THE MASS OF THE NEUTRAL \bar{K}
	mLc = 1.3828D0   !Mass of the Sigma*^+
!!*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
	g_DLcN = DSQRT (0.0372D0 * 4*pi)    ! as N-Y*-K from Holzenkamp et al.
	gg = 2.D0* g_DLcN**2 / (16.D0 * pi**2)  !2 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          ENDIF

	IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
	  fypiD = 0.D0
	  RETURN
	ENDIF

	ss0 = 0.D0
	kTmax = 10.D0
	kTint = kTmax/1000.D0

        DO ikT = 0, 1000
            kT = kTint * DBLE(ikT)

          kT2 = kT**2
	  SLcD = (kT2 + mD0**2)/(1.D0-y) + (kT2 + mLc**2)/y
!********************************************************************************
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sLcD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = DEXP( (mN**2 - sLcD)/(2*L**2) )        ! expon -- HSS
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - sLcD)/(L**2) )          ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mD0**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN                      ! DIPOLE -- s-channel Lambda exchange
            sM = (kT2 + (1.D0+y)*mD0**2.D0)/y + (kT2 
     &                         + y*mLc**2.D0)/(1.D0-y) + mN**2.D0
            FF = (L**4 + mLc**4)/(L**4 + sM**2) 
          ENDIF

      ss = ( kT2 + (mLc-y*mN)**2) * ( kT2 + (mLc+y*mN)**2 )**2
     &       / ( 6.D0*mLc**2*y**3) / ( (1.D0-y)*(SLcD - mN**2) )**2 
     &       * FF**2 * (2*kT)

	  IF (ikT/2*2.NE.ikT) THEN
	    ss0 = ss0 + 4*ss
	  ELSE IF (ikT/2*2.EQ.ikT) THEN
	    ss0 = ss0 + 2*ss
	  ENDIF

	ENDDO
	fypiD =  gg/(mD0**2) * (1.D0-y) / y * (kTint/3) * ss0
	RETURN   
	END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
