      SUBROUTINE GEN(I)
      IMPLICIT None
	Integer N
	Double Precision WPS,PIE,M,MT,MS,P,K
      COMMON/KIN/ WPS,PIE,M(10),MT(10),MS(10),P(10,4),K(10,4),N
	Double Precision X,AALAM
      COMMON/XVAR/ X(10)
*KEEP,RGPARAS.
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE

*KEND.
      Double Precision PSQ,EI,PI,COST,COSTM,SINT,PHI
      Integer I
      IF(I.EQ.N-1) THEN
         MT(I+1) = M(N)
      ELSE
         MT(I+1) = MS(I+1) + (MT(1)-MS(1))*X(3*I)
      ENDIF
      PSQ = AALAM(MT(I)**2,M(I)**2,MT(I+1)**2)/(4.D0*MT(I)**2)
      IF(PSQ.LT.0.0D0) WRITE(6,*) 'PSQ.LT.0 ',PSQ
      IF(PSQ.LT.0.0D0) PSQ=0.D0
      EI = DSQRT(PSQ+M(I)**2)
      PI = DSQRT(PSQ)
      COST = 1.D0-2.D0*X(3*I-1)
      IF(IPRO.EQ.21.OR.IPRO.EQ.20) THEN
      ELSE
c new check limits on PT2CUT --> make limits on cost
         if(pi**2.lt.pt2cut(ipro)) then
c         write(6,*) 'gen pi**2,pt2cut,ipro',pi**2,pt2cut(ipro),ipro
            wps = 0.d0
            return
         endif
         costm = dsqrt(1.d0 - pt2cut(ipro)/pi**2)
         cost = costm - 2.d0*X(3*I-1)*costm
c      write(6,*) 'gen: pi,ipro',pi,ipro
c      write(6,*) 'gen pt2cut,costm,cost',pt2cut(ipro),costm,cost,ipro
c new weight because of restricted phase space
         wps = wps*costm
      ENDIF
c end new
      SINT = DSQRT(1.D0-COST**2)
c new try for 1/pt**2 weighting (1/sin**2 )
cnew      SINT2 = X(3*I-1)
cnew      SINT = DSQRT(SINT2)
cnew      COST = DSQRT(1.D0 - SINT2)
cnew      IF(draprn().GT.0.5) COST = - COST
c      write(6,*) ' gen sin cos ',sint,cost
c end new
      PHI = 2.D0*PIE*X(3*I-2)
      P(I,1) = EI
      P(I,2) = PI*SINT*DSIN(PHI)
      P(I,3) = PI*SINT*DCOS(PHI)
c      write(6,*) ' GEN: PT**2 ',P(I,2)**2+P(I,3)**2,I
      P(I,4) = PI*COST
      if(p(I,2).ne.p(i,2)) then
         write(6,*) '  error in gen '
         write(6,*) ' pi,sint,phi ',pi,sint,phi,X(3*I-2),X(3*I-1)
      endif
      WPS = WPS*PI
      RETURN
      END
