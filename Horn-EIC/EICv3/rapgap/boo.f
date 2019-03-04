      SUBROUTINE BOO(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K	
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Integer I,J
	Double Precision BG1,G1,BG2,G2,BG3,G3,E,X,Z,Y
      BG1 = K(I,2)/MT(I)
      G1 = DSQRT(1.D0 + BG1**2)
      BG2 = K(I,3)/(MT(I)*G1)
      G2 = DSQRT(1.D0 + BG2**2)
      BG3 = K(I,4)/(MT(I)*G1*G2)
      G3 = DSQRT(1.D0 + BG3**2)
      E = G1*G2*G3*P(I,1) + BG1*G2*G3*P(I,2)
     +    + BG2*G3*P(I,3) + BG3*P(I,4)
      X = BG1*P(I,1) +G1*P(I,2)
      Y = BG2*G1*P(I,1) + BG2*BG1*P(I,2) +G2*P(I,3)
      Z = BG3*G2*G1*P(I,1) + BG3*G2*BG1*P(I,2)
     +    + BG3*BG2*P(I,3) +G3*P(I,4)
      P(I,1) = E
      P(I,2) = X
      P(I,3) = Y
      P(I,4) = Z
      IF(I.EQ.N-1) GOTO 20
      DO 10 J=1,4
   10 K(I+1,J) = K(I,J) - P(I,J)
C     WRITE(6,*) 'BOO ',(P(I,IK),IK=1,4)
      RETURN
   20 DO 30 J =1,4
   30 P(N,J) = K(I,J) - P(I,J)
      RETURN
      END
