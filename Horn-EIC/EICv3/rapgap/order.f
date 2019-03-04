      SUBROUTINE ORDER
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      IF(N.LE.3) RETURN
      WRITE(6,*) ' N .GT. 3; N = ',N
      STOP
      END
