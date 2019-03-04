	

      SUBROUTINE draprnv(RVEC,LENV)

      DOUBLE PRECISION TWOM24
      PARAMETER (TWOM24 = 2.D0**(-24))
	Integer lenv
	Double Precision RVEC(LENV)
	Real R
	Integer ld,li
	Parameter (ld=100)
	DIMENSION R(ld)
	if(lenv.gt.ld/2) then
	  write(6,*) ' draprnv : len > ld/2 ',len,ld/2
	  stop
	  endif
*  get two random numbers with single precision
      li=lenv*2
      CALL H1RNV (R, li)
      DO I=1,LENV
      RVEC(I) = DBLE(R(I)) + DBLE(R(LEN+I))*TWOM24
	enddo

	RETURN
      END
