*CMZ :  2.08/00 06/06/99  17.07.14  by  Hannes Jung
*CMZ :  2.00/01 20/04/95  17.18.33  by  Hannes Jung
*-- Author :
C*********************************************************************

      FUNCTION PYGAMM(X)

C...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
C...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
C...(Dover, 1965) 6.1.36.
      Double Precision B
      DIMENSION B(8)
      DATA B/-0.577191652D0,0.988205891D0,-0.897056937D0,0.918206857D0,
     +-0.756704078D0,0.482199394D0,-0.193527818D0,0.035868343D0/

      NX=INT(X)
      DX=X-NX

      PYGAMM=1.
      DXP=1.
      DO 10  I=1,8
         DXP=DXP*DX
   10 PYGAMM=PYGAMM+SNGL(B(I))*DXP
      IF(X.LT.1.) THEN
         PYGAMM=PYGAMM/X
      ELSE
         DO 20 IX=1,NX-1
   20    PYGAMM=(X-IX)*PYGAMM
      ENDIF
      RETURN
      END
