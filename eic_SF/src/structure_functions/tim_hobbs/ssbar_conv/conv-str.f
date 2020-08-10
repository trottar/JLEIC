C **********************************************************************
	PROGRAM conv_str
C
C  Computes the convolution of two distribution functions for the sake of
C  determining the final quark density of strange in the proton; the separate
C  distribution functions (the p --> L K vertex, and s-sbar distributions
C  in the Lambda or K are called by the leading program "convolute," and appear
C  as FORTRAN functions at the end of this file.
C
C  The choice of form factor, cutoffs, and constituent masses may be made
C  separately in the distribution function codes.
C
C  Written: T. Hobbs (May 10, 2012)
C  Modified: T. Hobbs (April 03, 2013)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,iy,iw,nw,FLAG,typ,dis
	PARAMETER (nx=95)
	PARAMETER (nw=95)
	EXTERNAL fyLcD, fycDb, f_DstL, f_DstSig, fypiD
	EXTERNAL fycLam, fycDvec, fycRS
	REAL*8  fycLam,fycLam0,fycLam1,fycDvec,fycDvec0,fypiD,LamQ,LamQB
	REAL*8  fyLcD, fycDb, f_DstL, f_DstSig, fycRS, fycRS0, fycRS1
	REAL*8  x,xmin,xmax,xint,y,ymax,y0,yint, intc
	REAL*8  w,wmin,wmax,wint,L,QDF_cL(nx),QDF_cS(nx),QDF_cbD(nx)
	REAL*8  xpt(nx),c(nx),c0,vDL(nx),c_2body(nx)
	REAL*8  BHPS(nx),vDLN0,vDLN(nw),QDF_cbDs(nx)
	REAL*8  fycDbN0,fycDbN1,fycDbFN(nw),fycDbN1,fycDbFN1(nw)
	REAL*8  fycLamN(nw),fycDvecN(nw),fycRSN(nw)
	REAL*8  fycLamN1(nw),fycRSN1(nw),HSF_RS1(nx),HSF_RS2(nx)
	REAL*8  HSF_PS1(nx),HSF_PS2(nx),HSF_PS3(nx),HSF_V1(nx),HSF_V2(nx)
C***********************************************************************
!---------------------------------------------------------------------
!WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
!      L = 2.98D0
      L = 1.D0
      LamQ = 1.D0
      LamQB = 8.D0
C***********************************************************************
!HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
!      typ = 1 !DIPOLE FORM FACTOR
      typ = 2 !EXPONENTIAL FORM FACTOR
C***********************************************************************
!THIS FLAG ALLOWS ONE TO SET A UNVERSAL CHOICE OF DISSOCIATION MODE
      dis = 0 
!---------------------------------------------------------------------
C***********************************************************************
!THIS CHOOSES AMONG THE VARIOUS POSSIBLE COMBINATIONS OF SPLITTING FUNCTION/PDF
      FLAG = 24
!THE FLAGS TOGGLE AMONG DISSOCIATION MODES AS: 
!     FLAG = 0,1  --- strange/ANTI-strange in Dbar^0 + Lambda^+_c  | J = 0 + 1/2
!     FLAG = 2,3  --- strange/ANTI-strange in Dbar^0 + Sigma^+_c   |
!     FLAG = 4,5  --- strange/ANTI-strange in Dbar^- + Sigma^++_c  |_________________
!     FLAG = 6,7  --- strange/ANTI-strange in Dbar^{*0} + Lambda^+_c  | J = 1 + 1/2
!     FLAG = 8,9  --- strange/ANTI-strange in Dbar^{*0} + Sigma^+_c   |
!     FLAG = 22/23 --- strange/ANTI-strange in D^{*-} + Sigma^++_c    |______________
!     FLAG = 10,11  --- strange/ANTI-strange in Dbar^{*0} + Sigma^{*+}_c  | J = 1 + 3/2
!     FLAG = 12,13  --- strange/ANTI-strange in Dbar^{-*} + Sigma^{*++}_c |____________
!     FLAG = 14,15  --- strange/ANTI-strange in Dbar^0 + Sigma^{*+}_c       | J = 0 + 3/2
!     FLAG = 16,17  --- strange/ANTI-strange in Dbar^- + Sigma^{*++}_c      |____________
!     FLAG = 18   --- incoherent sum: 2x * (c + sbar){x} = F2c      ****************** 
!     FLAG = 19   --- incoherent sum:  c(x)                         * F_2(x,Q=m_c) AND
!     FLAG = 20   --- incoherent sum:  sbar(x)                      * c-sbar_{TOT}
!     FLAG = 21   --- F2c as generated only with the dominant mode  ******************
!     FLAG = 24   --- x*{s-sbar}(x) using the lightest mode
!---------------------------------------------------------------------
!WE MUST FIRST ENSURE THAT THE INPUTS TO THE CONVOLUTION ARE FULLY
!NORMALIZED -- WE INTEGRATE THE EXTERNAL FUNCTIONS OVER A FULL
!RANGE 0 TO 1:
	iw = 0
	wmin = 0.01D0
	wmax = 1.D0
	wint = (wmax-wmin)/nw

       vDLN0 = 0.D0
       fycDbN0 = 0.D0
       fycDbN1 = 0.D0
       fycLam0 = 0.D0
       fycLam1 = 0.D0
       fycDvec0 = 0.D0
       fycRS0 = 0.D0
       fycRS1 = 0.D0

	DO   w = wmin, wmax, wint

	     vDLN(iw) = fyLcD(1.D0-w,L,typ,dis)
!	     vDLN(iw) = w*fyLcD(w,L,typ,dis)   !FOR MOMENTUM FRACTIONS

	     fycDbFN(iw) = fycDb(w,L,typ,0)
	     fycDbFN1(iw) = fycDb(w,L,typ,1)
	     fycLamN(iw) = fycLam(w,L,typ,0)
!	     fycLamN(iw) = fycLam(w,LamQ,typ,0)
	     fycLamN1(iw) = fycLam(w,L,typ,1)
	     fycDvecN(iw) = fycDvec(w,L,typ,0)
!	     fycDvecN(iw) = fycDvec(w,LamQB,typ,0)
	     fycRSN(iw) = fycRS(w,L,typ,0)
	     fycRSN1(iw) = fycRS(w,L,typ,1)

	  IF (iw/2*2.NE.iw) THEN
	     vDLN0 = vDLN0 + 4*vDLN(iw)
	     fycDbN0 = fycDbN0 + 4*fycDbFN(iw)
	     fycDbN1 = fycDbN1 + 4*fycDbFN1(iw)
	     fycLam0 = fycLam0 + 4*fycLamN(iw)
	     fycLam1 = fycLam1 + 4*fycLamN1(iw)
	     fycDvec0 = fycDvec0 + 4*fycDvecN(iw)
	     fycRS0 = fycRS0 + 4*fycRSN(iw)
	     fycRS1 = fycRS1 + 4*fycRSN1(iw)
	  ELSE IF (iw/2*2.EQ.iw) THEN
	     vDLN0 = vDLN0 + 2*vDLN(iw)
	     fycDbN0 = fycDbN0 + 2*fycDbFN(iw)
	     fycDbN1 = fycDbN1 + 2*fycDbFN1(iw)
	     fycLam0 = fycLam0 + 2*fycLamN(iw)
	     fycLam1 = fycLam1 + 2*fycLamN1(iw)
	     fycDvec0 = fycDvec0 + 2*fycDvecN(iw)
	     fycRS0 = fycRS0 + 2*fycRSN(iw)
	     fycRS1 = fycRS1 + 2*fycRSN1(iw)
	  ENDIF

	  iw = iw + 1

        ENDDO
	vDLN0 = (wint/3) * vDLN0
	fycDbN0 = (wint/3) * fycDbN0
	fycDbN1 = (wint/3) * fycDbN1
	fycLam0 = (wint/3) * fycLam0
	fycLam1 = (wint/3) * fycLam1
	fycDvec0 = (wint/3) * fycDvec0
	fycRS0 = (wint/3) * fycRS0
	fycRS1 = (wint/3) * fycRS1
!---------------------------------------------------------------------

! WE NOW COMPUTE CONVOLUTIONS AND DISTRIBUTIONS OVER A FULL RANGE OF x
	ix = 1
	xmin = 0.01D0
	xmax = 1.D0
	xint = xmax/nx
	c0 = 0.D0
	DO   x=xmin,xmax+xint/10,xint

	     xpt(ix) = x

! DEFINE THE BOUNDS OF THE CONVOLUTION INTEGRAL -- x to 1
	     c(ix) = 0.D0
	  iy = 0
	  y0 = x
	     ymax = 1.D0
	     yint = (ymax-y0)/1000
	  DO y=y0,ymax,yint

!!!!!!!!!!!!!! THE INTEGRAND OF THE CONVOLUTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!COMPUTE THE CONVOLUTION OF A HADRONIC VERTEX WITH A PARTON DISTRIBUTION FUNCTION

          IF (FLAG.EQ.0) THEN
            intc = (1 / y) * fyLcD(y,L,typ,0) * fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !s in Lambda_c + Dbar^0
          ELSE IF (FLAG.EQ.1) THEN
            intc = (1 / y) * fyLcD(1.D0-y,L,typ,0)*fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Lambda_c + Dbar^0
          ELSE IF (FLAG.EQ.2) THEN
            intc = (1 / y) * fyLcD(y,L,typ,1) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !s in Sigma_c + Dbar^0
          ELSE IF (FLAG.EQ.3) THEN
            intc = (1 / y) * fyLcD(1.D0-y,L,typ,1) * fycDb(x/y,L,typ,1)
     &                                         /fycDbN0 !sbar in Sigma_c + Dbar^0
          ELSE IF (FLAG.EQ.4) THEN
            intc = (1 / y) * fyLcD(y,L,typ,2) * fycLam(x/y,L,typ,2)
     &                                         /fycLam1 !s in Sigma^++_c + Dbar^-
          ELSE IF (FLAG.EQ.5) THEN
            intc = (1 / y) * fyLcD(1.D0-y,L,typ,2) * fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma^++_c + Dbar^-
          ELSE IF (FLAG.EQ.6) THEN
            intc = (1 / y)*f_DstL(1.D0-y,L,typ,0)*fycLam(x/y,L,typ,0)
!            intc = (1 / y)*f_DstL(1.D0-y,L,typ,0)*fycLam(x/y,LamQ,typ,0)
     &                                         /fycLam0 !s in Lambda_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.7) THEN
            intc = (1 / y) * f_DstL(y,L,typ,0) * fycDvec(x/y,L,typ,0)
!            intc = (1 / y) * f_DstL(y,L,typ,0) *fycDvec(x/y,LamQB,typ,0)
     &                                         /fycDvec0 !sbar in Lambda_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.8) THEN
            intc = (1 / y) * f_DstL(1.D0-y,L,typ,1)*fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !s in Sigma_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.9) THEN
            intc = (1 / y) * f_DstL(y,L,typ,1) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.10) THEN
            intc = (1 / y)*f_DstSig(1.D0-y,L,typ,0)*fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !s in Sigma*_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.11) THEN
            intc = (1 / y)*f_DstSig(y,L,typ,0)*fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Sigma*_c + Dbar^{*0}
          ELSE IF (FLAG.EQ.12) THEN
            intc = (1 / y)*f_DstSig(1.D0-y,L,typ,1)*fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !s in Sigma*_++c + D^{*-}
          ELSE IF (FLAG.EQ.13) THEN
            intc = (1 / y)*f_DstSig(y,L,typ,1)*fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma*++_c + D^{*-}
          ELSE IF (FLAG.EQ.14) THEN
            intc = (1 / y)*fypiD(y,L,typ,0)*fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !s in Sigma*+_c + Dbar^0
          ELSE IF (FLAG.EQ.15) THEN
            intc = (1 / y)*fypiD(1.D0-y,L,typ,0)*fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Sigma*+_c + Dbar^0
          ELSE IF (FLAG.EQ.16) THEN
            intc = (1 / y)*fypiD(y,L,typ,1)*fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !s in Sigma*++_c + D^-
          ELSE IF (FLAG.EQ.17) THEN
            intc = (1 / y)*fypiD(1.D0-y,L,typ,1)*fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma*++_c + D^-
          ELSE IF (FLAG.EQ.18) THEN
        intc = (4.D0/9.D0)*x*(1/y)*(fyLcD(y,L,typ,0)*fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !s in Lambda_c + Dbar^0
     &        + fyLcD(1.D0-y,L,typ,0) * fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Lambda_c + Dbar^0
     &        + fyLcD(y,L,typ,1) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !s in Sigma_c + Dbar^0
     &        + fyLcD(1.D0-y,L,typ,1) * fycDb(x/y,L,typ,1)
     &                                         /fycDbN0 !sbar in Sigma_c + Dbar^0
     &        + fyLcD(y,L,typ,2) * fycLam(x/y,L,typ,2)
     &                                         /fycLam1 !s in Sigma^++_c + Dbar^-
     &        + fyLcD(1.D0-y,L,typ,2) * fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma^++_c + Dbar^-
     &        + f_DstL(1.D0-y,L,typ,0) * fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !s in Lambda_c + Dbar^{*0}
     &        + f_DstL(y,L,typ,0) * fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Lambda_c + Dbar^{*0}
     &        + f_DstL(1.D0-y,L,typ,1) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !s in Sigma_c + Dbar^{*0}
     &        + f_DstL(y,L,typ,1) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c + Dbar^{*0}
     &        + f_DstSig(1.D0-y,L,typ,0) * fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !s in Sigma*_c + Dbar^{*0}
     &        + f_DstSig(y,L,typ,0)*fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Sigma*_c + Dbar^{*0}
     &        + f_DstSig(1.D0-y,L,typ,1) * fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !s in Sigma*_++c + Dbar^{*-}
     &        + f_DstSig(y,L,typ,1)*fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma*++_c + Dbar^{*-}
     &        + fypiD(y,L,typ,0)*fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !s in Sigma*+_c + Dbar^0
     &        + fypiD(1.D0-y,L,typ,0)*fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Sigma*+_c + Dbar^0
     &        + fypiD(y,L,typ,1)*fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !s in Sigma*++_c + D^-
     &        + fypiD(1.D0-y,L,typ,1)*fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma*++_c + D^-
     &        + f_DstL(1.D0-y,L,typ,2) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !s in Sigma_c^++ + D^{*-}
     &        + f_DstL(y,L,typ,2) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c^++ + D^{*-}
     &           )  !THE F2c STRUCTURE FUNCTION!
          ELSE IF (FLAG.EQ.19) THEN
            intc = (1 / y)*(fyLcD(y,L,typ,0) * fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !c in Lambda_c + Dbar^0
     &        + fyLcD(y,L,typ,1) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !c in Sigma_c + Dbar^0
     &        + fyLcD(y,L,typ,2) * fycLam(x/y,L,typ,2)
     &                                         /fycLam1 !c in Sigma^++_c + Dbar^-
     &        + f_DstL(1.D0-y,L,typ,0) * fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !c in Lambda_c + Dbar^{*0}
     &        + f_DstL(1.D0-y,L,typ,1) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !c in Sigma_c + Dbar^{*0}
     &        + f_DstSig(1.D0-y,L,typ,0) * fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !c in Sigma*_c + Dbar^{*0}
     &        + f_DstSig(1.D0-y,L,typ,1) * fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !c in Sigma*_++c + Dbar^{*-}
     &        + fypiD(y,L,typ,0)*fycRS(x/y,L,typ,0)
     &                                         /fycRS0 !c in Sigma*+_c + Dbar^0
     &        + fypiD(y,L,typ,1)*fycRS(x/y,L,typ,1)
     &                                         /fycRS1 !c in Sigma*++_c + D^-
     &        + f_DstL(1.D0-y,L,typ,2) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1 !c in Sigma_c^++ + D^{*-}
     &           )  !THE INCOHERENT SUM -- CHARM IN THE PROTON

          ELSE IF (FLAG.EQ.20) THEN
           intc = (1 / y)*(fyLcD(1.D0-y,L,typ,0) * fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Lambda_c + Dbar^0
     &        + fyLcD(1.D0-y,L,typ,1) * fycDb(x/y,L,typ,1)
     &                                         /fycDbN0 !sbar in Sigma_c + Dbar^0
     &        + fyLcD(1.D0-y,L,typ,2) * fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma^++_c + Dbar^-
     &        + f_DstL(y,L,typ,0) * fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Lambda_c + Dbar^{*0}
     &        + f_DstL(y,L,typ,1) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c + Dbar^{*0}
     &        + f_DstSig(y,L,typ,0)*fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Sigma*_c + Dbar^{*0}
     &        + f_DstSig(y,L,typ,1)*fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma*++_c + Dbar^{*-}
     &        + fypiD(1.D0-y,L,typ,0)*fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 !sbar in Sigma*+_c + Dbar^0
     &        + fypiD(1.D0-y,L,typ,1)*fycDb(x/y,L,typ,2)
     &                                         /fycDbN1 !sbar in Sigma*++_c + D^-
     &        + f_DstL(y,L,typ,2) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c^++ + D^{*-}
     &           )  !THE INCOHERENT SUM -- ANTI-CHARM IN THE PROTON
          ELSE IF (FLAG.EQ.21) THEN
        intc = (4.D0/9.D0)*x*(1/y)*(
     &        + f_DstL(1.D0-y,L,typ,0) * fycLam(x/y,L,typ,0)
     &                                         /fycLam0 !c in Lambda_c + Dbar^{*0}
     &        + f_DstL(y,L,typ,0) * fycDvec(x/y,L,typ,0)
     &                                         /fycDvec0 !sbar in Lambda_c + Dbar^{*0}
     &           )  !THE F2c USING ONLY THE DOMINANT MODE
          ELSE IF (FLAG.EQ.22) THEN
        intc = (1/y) * f_DstL(1.D0-y,L,typ,2) * fycLam(x/y,L,typ,1)
     &                                         /fycLam1  !c in Sigma_c^++ + D^{*-}
          ELSE IF (FLAG.EQ.23) THEN
        intc = (1/y) * f_DstL(y,L,typ,2) * fycDvec(x/y,L,typ,1)
     &                                         /fycDvec0 !sbar in Sigma_c^++ + D^{*-}
          ELSE IF (FLAG.EQ.24) THEN
            intc = 4.5D0*14.423D0 *x * ( (1 / y)*fyLcD(y,L,typ,0) 
     &                   * fycLam(x/y,L,typ,0)/fycLam0
     &           - (1 / y)*fyLcD(1.D0-y,L,typ,0)*fycDb(x/y,L,typ,0)
     &                                         /fycDbN0 ) ! x*{s-sbar}(x)
          ENDIF
        
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
	     IF (iy.EQ.0) THEN
	      c(ix) = c(ix) + intc
	    ELSE IF (iy/2*2.NE.iy) THEN
	      c(ix) = c(ix) + 4*intc
	    ELSE IF (iy/2*2.EQ.iy) THEN
	      c(ix) = c(ix) + 2*intc
	     ENDIF

 123	    iy = iy + 1

	  ENDDO
! AS A CHECK, WE ALSO WILL PRINT THE INPUT DISTRIBUTIONS TO A SEPARATE FILE:
!RECIP
!	     HSF_PS1(ix) = fyLcD(1.D0-x,L,typ,0)
	     HSF_PS1(ix) = fyLcD(x,L,typ,0)
	     HSF_PS2(ix) = fyLcD(1.D0-x,L,typ,1)
	     HSF_PS3(ix) = fyLcD(1.D0-x,L,typ,2)
	     HSF_V1(ix) = f_DstL(x,L,typ,0)
	     HSF_V2(ix) = f_DstL(x,L,typ,1)
	     HSF_RS1(ix) = f_DstSig(x,L,typ,0)
	     HSF_RS2(ix) = f_DstSig(x,L,typ,1)

	     QDF_cL(ix) = fycLam(x,L,typ,0)/fycLam0
	     QDF_cS(ix) = fycRS(x,L,typ,0)/fycRS0
	     QDF_cbD(ix) = fycDb(x,L,typ,0)/fycDbN0
	     QDF_cbDs(ix) = fycDvec(x,L,typ,0)/fycDvec0

! THE BHPS RESULT IS COMPUTED HERE ALSO, FOR COMPARISON
	     BHPS(ix) = 18.D0 * x**2 * ( (1./3.)*(1.D0 - x)*(1.D0 + 10.D0*x
     &           + x**2) - 2.D0 * x * (1.D0 + x) *DLOG(1.D0/x) )

	     c(ix) = (yint/3) * c(ix)

	     c(1) = 0.D0

	print *, 'ix,x,cbD,cbD*,cS*=',ix,x,QDF_cbD(ix),QDF_cbDs(ix),QDF_cS(ix)

	  IF (ix/2*2.NE.ix) THEN
	     c0 = c0 + 4*c(ix)
	  ELSE IF (ix/2*2.EQ.ix) THEN
	     c0 = c0 + 2*c(ix)
	  ENDIF

	  ix = ix + 1

	ENDDO
	c0 = (xint/3) * c0

	print *, 'First moment of sbar(x)=',c0
	print *, 'First moment of physical D_{M/B}=', vDLN0
	print *, 'First moment of spinor QM sbar in D^0=', fycDbN0

c       if(.true.)stop
C...Write to file
          IF (FLAG.EQ.0) THEN
	OPEN (11,FILE='DAT-str/strng_LcD_EM.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.1) THEN
	OPEN (11,FILE='DAT-str/sbar_LcD_EM.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.2) THEN
	OPEN (11,FILE='DAT-str/strange_ScD.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.3) THEN
	OPEN (11,FILE='DAT-str/sbar_ScD.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.4) THEN
	OPEN (11,FILE='DAT-str/strange_Sc++D-.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.5) THEN
	OPEN (11,FILE='DAT-str/sbar_Sc++D-.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.6) THEN
	OPEN (11,FILE='DAT-str/strng_LKs_EM.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.7) THEN
	OPEN (11,FILE='DAT-str/sbar_LKs_EM.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.8) THEN
	OPEN (11,FILE='DAT-str/strange_ScDs.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.9) THEN
	OPEN (11,FILE='DAT-str/sbar_ScDs.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.10) THEN
	OPEN (11,FILE='DAT-str/strange_Ss+Ds.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.11) THEN
	OPEN (11,FILE='DAT-str/sbar_Ss+Ds.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.12) THEN
	OPEN (11,FILE='DAT-str/strange_Ss++Ds.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.13) THEN
	OPEN (11,FILE='DAT-str/sbar_Ss++Ds.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.14) THEN
	OPEN (11,FILE='DAT-str/strange_Ss+Db.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.15) THEN
	OPEN (11,FILE='DAT-str/sbar_Ss+Db.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.16) THEN
	OPEN (11,FILE='DAT-str/strange_Ss++Db.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.17) THEN
	OPEN (11,FILE='DAT-str/sbar_Ss++Db.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.18) THEN
	OPEN (11,FILE='DAT-str/F2s_tot-HQ.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.19) THEN
	OPEN (11,FILE='DAT-str/strange_sum_HQ.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.20) THEN
	OPEN (11,FILE='DAT-str/sbar_sum_HQ.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.21) THEN
	OPEN (11,FILE='DAT-str/F2s_dom_HQ.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.22) THEN
	OPEN (11,FILE='DAT-str/strange_S++Ds-.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.23) THEN
	OPEN (11,FILE='DAT-str/sbar_S++Ds-.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ELSE IF (FLAG.EQ.24) THEN
	OPEN (11,FILE='DAT-str/xs-sbar_EM-4p5%.dat',STATUS='UNKNOWN',
     &                                        FORM='FORMATTED')
          ENDIF
	  DO ix=1,nx
	  WRITE (11,*) xpt(ix), c(ix)
	ENDDO
	CLOSE (11)
	OPEN (12,FILE='Db_2p91_M2.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (12,*) xpt(ix), vDL(ix)
	ENDDO
	CLOSE (12)
	OPEN (13,FILE='inquark_dist.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (13,*) xpt(ix),c_2body(ix)
	ENDDO
	CLOSE (13)
	OPEN (14,FILE='BHPS.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (14,*) xpt(ix), (8.D0/9.D0)*xpt(ix)*(BHPS(ix))
	ENDDO
	CLOSE (14)
	OPEN (15,FILE='HSF_PS_LK.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (15,*) xpt(ix), HSF_PS1(ix)!, HSF_PS2(ix), HSF_PS3(ix)
	ENDDO
	CLOSE (15)
	OPEN (16,FILE='HSF_V.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (16,*) xpt(ix), HSF_V1(ix), HSF_V2(ix)
	ENDDO
	CLOSE (16)
	OPEN (17,FILE='HSF_RS.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (17,*) xpt(ix), HSF_RS1(ix), HSF_RS2(ix)
	ENDDO
	CLOSE (17)
	OPEN (18,FILE='QDF_s_EM.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (18,*) xpt(ix), QDF_cL(ix)!, QDF_cS(ix)
	ENDDO
	CLOSE (18)
	OPEN (19,FILE='QDF_sb_EM.dat',STATUS='UNKNOWN',FORM='FORMATTED')
	  DO ix=1,nx
	  WRITE (19,*) xpt(ix), QDF_cbD(ix)!, QDF_cbDs(ix)
	ENDDO
	CLOSE (19)
	END
C ***********************************************************************
	FUNCTION fyLcD (y,L,typ,dis)
C
C  Function giving numerical value of f(y) for N-Lambda-K vertex
C    y is l.c. momentum fraction on baryon.
C
C  Written: W. Melnitchouk (1999)
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
	ikT = 1
	kTmax = 10.D0
	kTint = kTmax/1000.D0
	DO kT=kTint,kTmax,kTint
	  kT2 = kT**2
	  SLcD = (kT2 + mD0**2)/(1.D0-y) + (kT2 + mLc**2)/y

          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + SLcD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + SLcD))**2    ! dipole
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
	  ikT = ikT + 1
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
	g_RoNN = DSQRT (0.841D0 * 4.D0 * pi * (2.D0)) !2 ISOSPIN
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
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          SRoN = (kT2 + mD**2)/y + (kT2 + mL**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))**2    ! dipole
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

!          ss = gg*sv                                         !VECTOR
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)

!          ss = fg*si/mN                                      !INTERFERENCE
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)

!          ss = fff*st/mN**2                                  !TENSOR
!     &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)

          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
          ikT = ikT + 1
        ENDDO
        f_DstL = y / (1.D0-y) * (kTint/3.D0) * ss0 
        RETURN
        END
C ***************************************************************************
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
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          
          SRoN = (kT2 + mD**2)/y + (kT2 + mS**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + sRoN))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mN**2 - SRoN)/(L**2) )        ! expon
          ELSE IF (typ.EQ.3) THEN
            t = ( -kT2 - (1.D0-y)*(mD**2 - y*mN**2) ) / y
            FF = ((L**2 - mS**2) / (L**2 - t))**2      ! t-channel dip
          ELSE IF (typ.EQ.4) THEN
	  sM = (kT2 + (2.D0-y)*mD**2.D0)/(1.D0-y) + (kT2 + (1.D0-y)*mS**2.D0)
     &       /y + mN**2.D0
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
          ikT = ikT + 1
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
	REAL*8  pi,mN,mD0,mLc,g_DLcN,gg,FF,L,sM

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
	ikT = 1
	kTmax = 10.D0
	kTint = kTmax/1000.D0
	DO kT=kTint,kTmax,kTint
	  kT2 = kT**2
	  SLcD = (kT2 + mD0**2)/(1.D0-y) + (kT2 + mLc**2)/y
!********************************************************************************
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mN**2) / (L**2 + sLcD))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mN**2) / (L**2 + sLcD))**2    ! dipole
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
	  ikT = ikT + 1
	ENDDO
	fypiD =  gg/(mD0**2) * (1.D0-y) / y * (kTint/3) * ss0
	RETURN   
	END
C ***************************************************************************
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C ***************************************************************************
C&&&&&&&&&  THE MODEL QUARK DISTRIBUTIONS BEGIN HERE  &&&&&&&&&&&&&&&&&&&&&&C
! ***************************************************************************
        FUNCTION fycLam (y,L,typ,dis)
!  Function giving numerical value of f(y) for the distribution of strange
!  inside the Lambda_c baryon
! ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SCoM,t,
     &          sv,ss,ss0
        REAL*8  y,fycLam
        REAL*8  pi,mN,mRo,FF,L
        REAL*8  mc,mL,mD
     
        pi = 4.D0 * DATAN(1.D0)

        mN = 0.93891897D0  !in GeV

        mc = 0.45D0     !LARGE strange quark mass

        mD = 1.D0  !Mass of the scalar (ud) diquark
          IF (dis.EQ.0) THEN
	mL = 1.1157D0  !THE MASS OF THE LAMBDA^0
          ELSE IF (dis.EQ.1) THEN
	mL = 1.1157D0  !THE MASS OF THE LAMBDA^0
          ELSE IF (dis.EQ.2) THEN
	mL = 1.1157D0  !THE MASS OF THE LAMBDA^0
          ENDIF


        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          fycLam = 0.D0
          RETURN
        ENDIF
     
        ss0 = 0.D0
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          
          SCoM = (kT2 + mc**2)/y + (kT2 + mD**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mL**2) / (L**2 + sCoM))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mL**2) / (L**2 + sCoM))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
!            FF = DEXP( (mL**2 - SCoM)/(L**2) )*(SCoM - mL**2)   ! expon
            FF = DEXP( (mL**2 - SCoM)/(L**2) )          ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mRo**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN
            IF (sCoM.LE.L**2) FF = 1.D0                 ! sharp cutoff
            IF (sCoM.GT.L**2) FF = 0.D0
          ENDIF
!***********************************************************************
!...TOPT vertex factor from the scalar-spinor-spinor coupling
          sv = ( kT2 + (y*mL + mc)**2.D0 ) / y
!***********************************************************************
          if (kt.gt.0.0 .and. (sv.lt.0.0) ) then
            print *,'CS3 -- ##### kT,y,sv,st =',kt,y,sv
            stop
          endif

          ss = sv
     &       / ( (SCoM - mL**2) )**2 * FF**2  * (2.D0*kT)

          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
          ikT = ikT + 1
        ENDDO

        fycLam = 1.D0 /y /(1.D0-y)*(kTint/3.D0) * ss0
        RETURN
        END

! ***************************************************************************
        FUNCTION fycDb (y,L,typ,dis)
!  Function giving numerical value of f(y) for the distribution of the sbar
!  inside the K meson
! ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SCoM,t,
     &          sv,ss,ss0
        REAL*8  y,fycDb
        REAL*8  pi,mN,mRo,FF,L
        REAL*8  mc,mu,mD
     
        pi = 4.D0 * DATAN(1.D0)

        mN = 0.93891897D0  !in GeV

        mc = 0.45D0     !LARGE anti-strange quark mass

        mu = mN / 3.D0
          IF (dis.EQ.0) THEN
        mD = 0.4937D0  !THE MASS OF THE K+
          ELSE IF (dis.EQ.1) THEN
        mD = 0.4937D0  !THE MASS OF THE K+
          ENDIF

        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          fycDb = 0.D0
          RETURN
        ENDIF
     
        ss0 = 0.D0
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          
          SCoM = (kT2 + mc**2)/y + (kT2 + mu**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mD**2) / (L**2 + sCoM))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mD**2) / (L**2 + sCoM))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
            FF = DEXP( (mD**2 - SCoM)/(L**2) )          ! expon
!            FF = DEXP( (mD**2 - SCoM)/(L**2) )*(SCoM - mD**2)  ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mRo**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN
            IF (sCoM.LE.L**2) FF = 1.D0                 ! sharp cutoff
            IF (sCoM.GT.L**2) FF = 0.D0
          ENDIF
!************************************************************************
!...TOPT vertex factor from the pseudoscalar-spinor-spinor coupling
          sv = ( kT2 + (y*mu + (1.D0 - y)*mc)**2.D0 )
     &         / (y*(1-y))
!************************************************************************
          if (kt.gt.0.0 .and. (sv.lt.0.0) ) then
            print *,'CS4 -- ##### kT,y,sv,st =',kt,y,sv
            stop
          endif
        
          ss = sv
     &       / ( (SCoM - mD**2) )**2 * FF**2  * (2.D0*kT)

          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
          ikT = ikT + 1
        ENDDO

        fycDb = 1.D0 /y /(1.D0-y)*(kTint/3.D0) * ss0
        RETURN
        END
! ***************************************************************************
! ***************************************************************************
        FUNCTION fycDvec (y,L,typ,dis)
!  Function giving numerical value of f(y) for the distribution of the c_bar
!  inside the spin-1 Dbar^{*0} strange meson
! ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SCoM,t
        REAL*8  y,fycDvec,sv,ss,ss0
        REAL*8  pi,mN,mRo,FF,L
        REAL*8  mc,mu,mD
     
        pi = 4.D0 * DATAN(1.D0)

        mN = 0.93891897D0  !proton mass in GeV
!*****************************************************************
!        mc = 1.75D0     !LARGE anti-strange quark mass
        mc = 0.45D0     !LARGE anti-strange quark mass
!*****************************************************************
        mu = mN / 3.D0   !Mass of the spectator quark
          IF (dis.EQ.0) THEN
        mD = 0.8917D0  !THE MASS OF THE K*+
          ELSE IF (dis.EQ.1) THEN
        mD = 0.8917D0  !THE MASS OF THE K*+
          ELSE IF (dis.EQ.2) THEN
        mD = 0.8917D0  !THE MASS OF THE K*+
          ENDIF
!*****************************************************************
!*****************************************************************
        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          fycDvec = 0.D0
          RETURN
        ENDIF
        ss0 = 0.D0
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          
          SCoM = (kT2 + mc**2)/y + (kT2 + mu**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mD**2) / (L**2 + sCoM))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mD**2) / (L**2 + sCoM))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
!            FF = DEXP( (mD**2 - SCoM)/(L**2) )*(SCoM - mD**2)  ! expon
            FF = DEXP( (mD**2 - SCoM)/(L**2) )         ! expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mRo**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN
            IF (sCoM.LE.L**2) FF = 1.D0                 ! sharp cutoff
            IF (sCoM.GT.L**2) FF = 0.D0
          ENDIF
C********************************************************************
!...TOPT vertex factor from the (pseudo)vector-spinor-spinor coupling
!... ASSUMING A VECTOR D MESON .................................!
          sv = (  ((mu**2+kT2)/(mD**2) + (1.D0-y)**2)*(mc**2+kT2
     &       + y**2*mD**2) + kT2 + y**2*mu**2 + (1.D0-y)**2*mc**2
     &       + 6.D0*y*(1.D0-y)*mu*mc   )
     &         / (y*(1-y))
C********************************************************************
C********************************************************************
          if (kt.gt.0.0 .and. (sv.lt.0.0) ) then
            print *,'CS5 -- ##### kT,y,sv,st =',kt,y,sv
            stop
          endif
        
          ss = sv
     &       / ( (SCoM - mD**2) )**2 * FF**2  * (2.D0*kT)

          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
          ikT = ikT + 1
        ENDDO

        fycDvec = 1.D0 /y /(1.D0-y)*(kTint/3.D0) * ss0
        RETURN
        END
! ***************************************************************************
! ***************************************************************************
        FUNCTION fycRS (y,L,typ,dis)
!  Function giving numerical value of f(y) for the distribution of the strange quark
!  inside the spin-3/2 Sigma* strange baryon
! ***************************************************************************
        IMPLICIT NONE
        INTEGER ikT,typ,dis
        REAL*8  kT,kT2,kTmax,kTint,SCoM,t,sv,ss,ss0
        REAL*8  y,fycRS,pi,mN,mRo,FF,L
        REAL*8  mc,mS,mD,P_p,P_k,k_p,De,PD,kD,p_D

        pi = 4.D0 * DATAN(1.D0)
!---------------------------------------------------------------------
!**************  MASS CHOICES ****************************************
!---------------------------------------------------------------------
        mN = 0.93891897D0  !proton mass in GeV
!********************************************************************
        mc = 0.45D0         !LARGE anti-strange quark mass
!********************************************************************
        mD = 1.D0              !Mass of the axial-vector diqaurk
          IF (dis.EQ.0) THEN
        mS = 2.5175D0          !Mass of the strange Sigma*^+
          ELSE IF (dis.EQ.1) THEN
        mS = 2.5184D0          !Mass of the strange Sigma*^++
          ENDIF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
          fycRS = 0.D0
          RETURN
        ENDIF
        ss0 = 0.D0
        ikT = 0
        kTmax = 10.D0
        kTint = kTmax/1000.D0  
        DO kT=0.D0,kTmax,kTint
          kT2 = kT**2
          
          SCoM = (kT2 + mc**2)/y + (kT2 + mD**2)/(1.D0-y)
        
          IF (typ.EQ.0) THEN
            FF = ((L**2 + mS**2) / (L**2 + sCoM))       ! monopole
          ELSE IF (typ.EQ.1) THEN
            FF = ((L**2 + mS**2) / (L**2 + sCoM))**2    ! dipole
          ELSE IF (typ.EQ.2) THEN
!            FF = DEXP( (mS**2 - SCoM)/(L**2) ) * (SCoM - mS**2) !expon
            FF = DEXP( (mS**2 - SCoM)/(L**2) )  !expon
          ELSE IF (typ.EQ.3) THEN
            t = (- kT2 - mN**2*y**2) /(1.D0-y)
            FF = ((L**2 - mRo**2) / (L**2 - t))**2      ! cov dip
          ELSE IF (typ.EQ.4) THEN
            IF (sCoM.LE.L**2) FF = 1.D0                 ! sharp cutoff
            IF (sCoM.GT.L**2) FF = 0.D0
          ENDIF

!...TOPT vertex factor from the RS-vector-spinor coupling...
!-----------------------------------------------------------------------
!...THE LORENTZ-INVARIANT INNER PRODUCTS ARE DEFINED:
          P_p = (kT2 + mD**2 + (1.D0-y)**2*mS**2)/(2*(1.D0-y))
          P_k = (kT2 + mc**2 + y**2*mS**2)/(2*y)
          k_p = (kT2 + y**2*mD**2 + (1.D0-y)**2
     &                                     *mc**2)/(2*y*(1.D0-y))

          De = mS**2.D0 + mc**2 - 2.D0*P_k  !='DELTA'**2
          PD = mS**2.D0 - P_k
          kD = P_k - mc**2
          p_D = P_p - k_p

!...FOLLOWED BY THE EXTENDED VERTEX FACTOR COMPUTED FROM IMF DIAGRAMS
!NO SIMPLIFICATION --- CODED DIRECTLY FROM NOTES ----------
!------ THE PSEUDOVECTOR TRACE FACTORS ----------
         sv = (-16*(De*P_k - (3./2.)*De*mc*mS - PD*kD)   !TERMS I+II
     &  + (8./mD**2)*(-p_D*(p_D*P_k - PD*k_p - kD*P_p) - De*k_p*P_p 
     &  + mc*mS*p_D**2) + 8.*(-kD*PD + De*P_k - mc*mS*De))/3.D0

     &  - (1./(3*mS))*(-8*mc*(mS**2*De + PD**2)   !TERM III
     &  + (8*mc/(mD**2))*((-2*PD*p_D*P_p) + (De*P_p**2)))

     &  - (2./(3.*mS**2))*(-8.*PD**2*P_k   !TERM IV
     &  - 4*De*mS**2*P_k + 4.*mS**3*mc*De + 8*mc*mS*PD**2) 
     &  - (2./(3.*mD**2*mS**2))*(4*mD**2*PD**2-8*PD*p_D*P_p
     &  + 4*De*P_p**2)*(P_k-mc*mS)
!-----------------------------------------------------------------------
      IF (kt.gt.0.0 .AND. (sv.lt.0.0)) THEN
        print *,'CS6 -- ##### kT,y,sv =',kt,y,sv
        STOP
      ENDIF
      ss = (1./(6.*pi**2*mD**2)) * sv      !FROM NOTES
     &       / ( (SCoM - mS**2) )**2 * FF**2  * (2.D0*kT)

          IF (ikT.EQ.0) THEN
            ss0 = ss0 + ss
          ELSE IF (ikT/2*2.NE.ikT) THEN
            ss0 = ss0 + 4.D0*ss
          ELSE IF (ikT/2*2.EQ.ikT) THEN
            ss0 = ss0 + 2.D0*ss
          ENDIF
          
          ikT = ikT + 1
        ENDDO

        fycRS = 1.D0 /y /(1.D0-y)*(kTint/3.D0) * ss0
        RETURN
        END
! ***************************************************************************
! ***************************************************************************
! ***************************************************************************
