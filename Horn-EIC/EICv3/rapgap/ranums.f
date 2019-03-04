      SUBROUTINE RANUMS (X,N)
	IMPLICIT NONE
      Double Precision X(20)
	Double Precision VEC(20)
	Integer N,LEN,I
      LEN=N
      CALL draprnV(VEC,LEN)
      DO 10 I=1,N
   10 X(I)=VEC(I)
      RETURN
      END
