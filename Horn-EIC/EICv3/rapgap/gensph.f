*CMZ :  2.08/04 22/12/99  15.39.29  by  Hannes Jung
*CMZ :  2.06/34 20/05/98  09.12.55  by  Hannes Jung
*CMZ :  2.05/28 14/07/97  11.32.12  by  Hannes Jung
*CMZ :  2.01/09 10/02/96  18.32.43  by  Hannes Jung
*CMZ :  2.01/04 24/01/96  16.10.48  by  Hannes Jung
*CMZ :  2.01/03 20/01/96  17.13.34  by  Hannes Jung
*CMZ :  2.01/01 09/01/96  14.04.24  by  Hannes Jung
*-- Author :

      SUBROUTINE GENSPH(VEC)
C--------------------------------------------------------C
C GENERATES AN UNIT-LENGTH VECTOR UNIFORMLY IN 4PI.      C
C USES 1 SQRT AND NO SIN,COS'S. HITOSHI OCT-86           C
C--------------------------------------------------------C
      REAL VEC(3)
C
      Double Precision RNDMV(2)
   10 CALL draprnV(RNDMV,2)
      U1=2.*RNDMV(1)-1.
      U2=2.*RNDMV(2)-1.
      S=U1*U1+U2*U2
      IF(S.GE.1.) GOTO 10
C
      CO2=2.*SQRT(1.-S)
      VEC(1)=U1*CO2
      VEC(2)=U2*CO2
      VEC(3)=1.-2.*S
C
      RETURN
      END
