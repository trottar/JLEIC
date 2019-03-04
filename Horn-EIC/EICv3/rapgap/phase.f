*CMZ :  2.08/04 22/12/99  15.39.26  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.22.33  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  11.39.23  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      SUBROUTINE PHASE(NP,ET,XM,PCM,WT)
c to generate n-body phase space
c taken from: Collider Physics
c             V.D. Barger and R.J.N Phillips
c             Adison-Wesley Publishing Company, Inc.
c             1987
	Implicit none
	Integer NP
	Double Precision ET,XM,PCM,WT
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
      DIMENSION XM(18),PCM(4,18)
	Integer I,J,Jmax
      N=NP
      PIE=4.D0*DATAN(1.D0)
      WPS=0.D0
      DO 10 I=1,N
         M(I)=0.D0
         MT(I)=0.D0
         MS(I)=0.D0
         DO 10 J=1,4
            K(I,J)=0.D0
            P(I,J)=0.D0
   10 CONTINUE
      DO 20 I=1,N
         M(I)=XM(I)
   20 CONTINUE
      K(1,1)=ET
      MT(1) =ET
      DO 30 I=1,N
   30 MS(1) = MS(1) + M(I)
      CALL ORDER
      JMAX = N - 1
      WPS = JMAX*((MT(1)-MS(1))/(4.D0*PIE**2))**(N-2)/
     +      (4.D0*PIE*MT(1))
      DO 40 J =1,JMAX
         WPS = WPS/J
         MS(J+1) = MS(J) - M(J)
         CALL GEN(J)
   40 CALL BOO(J)
      DO 50 J=1,N
         PCM(1,J) = P(J,2)
         PCM(2,J) = P(J,3)
         PCM(3,J) = P(J,4)
         PCM(4,J) = P(J,1)
   50 CONTINUE
      if(p(1,2).ne.p(1,2)) write(6,*) ' error in phase'
      WT=WPS*(2.D0*PIE)**(3*N-4)
      RETURN
      END
