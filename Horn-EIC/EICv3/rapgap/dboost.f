*CMZ :  2.06/34 20/05/98  09.12.55  by  Hannes Jung
*CMZ :  2.01/09 10/02/96  18.32.43  by  Hannes Jung
*CMZ :  2.01/03 20/01/96  17.13.34  by  Hannes Jung
*CMZ :  2.01/01 09/01/96  14.04.24  by  Hannes Jung
*-- Author :

      SUBROUTINE DBOOST(P0,U0,PA,PB)
C---------------------------------------------------------C
C 4-MOMENTUM PA IN THE REST FRAME OF PARTICLE 0 IS        C
C EVALUATED IN THE FRAME WHERE PARTILE 0 HAS A 4-MOMENTUM C
C P0 AND RETURNED AS PB. U0 IS THE MASS OF THE PARTICLE 0.C
C PA AND PB CAN SHARE THE SAME MEMORY. OCT-86 HITOSHI     C
C THIS IS THE DOUBLE PRECISION VERSION OF LBOOST.         C
C---------------------------------------------------------C
      DOUBLE PRECISION P0(4),PA(4),PB(4),U0
      DOUBLE PRECISION EP,CONS
C
      EP=(PA(4)*P0(4)+PA(1)*P0(1)+PA(2)*P0(2)+PA(3)*P0(3))/U0
      CONS=(PA(4)+EP)/(U0+P0(4))
      PB(1)=PA(1)+CONS*P0(1)
      PB(2)=PA(2)+CONS*P0(2)
      PB(3)=PA(3)+CONS*P0(3)
      PB(4)=EP
C
      RETURN
      END
