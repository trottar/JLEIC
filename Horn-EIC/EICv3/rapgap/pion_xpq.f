*CMZ :  2.03/12 11/10/96  13.36.58  by  Hannes Jung
*CMZ :  2.03/08 18/09/96  19.03.02  by  Julian Phillips
*CMZ :  1.00/03 08/07/96  16.50.06  by  Julian Phillips
*-- Author :    Julian Phillips   08/07/96
C Construct PDFs for pion or other meson
      SUBROUTINE PION_XPQ(XX,Q2,XPQ)
      IMPLICIT REAL*8 (A-G,O-Z)
      DIMENSION XPQ(-6:6)
      DOUBLE PRECISION XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU

C PDFlib only contains data for the negative pion, so here the dval and
C uval are swapped to give a positive pion
*! ud isospin 0 meson
      ISF=2
C      ISF=1 ! pi+

      QQ = DSQRT(Q2)

C Initialise
      DO I=-6,6
         XPQ(I) = 0.D0
      ENDDO

      IF(XX.GE.1. .OR. XX.LE.0.) RETURN

      CALL STRUCTF(XX,QQ,UPV,DNV,SEA,STR,CHM,BOT,TOP,GLU)

C Positive Pion
      IF(ISF.EQ.1) THEN
*! Gluon
         XPQ(0)= GLU
*! D
         XPQ(1)= UPV+SEA
*! Dbar
         XPQ(-1)= SEA
*! U
         XPQ(2)= DNV+SEA
*! Ubar
         XPQ(-2)= SEA
*! S
         XPQ(3)= STR
*! Sbar
         XPQ(-3)= STR
*! C
         XPQ(4)= CHM
*! Cbar
         XPQ(-4)= CHM
*! B
         XPQ(5)= BOT
*! Bbar
         XPQ(-5)= BOT
*! T
         XPQ(6)= TOP
*! Tbar
         XPQ(-6)= TOP
      ELSEIF(ISF.EQ.2) THEN
C Flavour singlet meson (f2)
         VAL = UPV/2.D0
*! Gluon
         XPQ(0)= GLU
*! D
         XPQ(1)= VAL+SEA
*! Dbar
         XPQ(-1)= VAL+SEA
*! U
         XPQ(2)= VAL+SEA
*! Ubar
         XPQ(-2)= VAL+SEA
*! S
         XPQ(3)= STR
*! Sbar
         XPQ(-3)= STR
*! C
         XPQ(4)= CHM
*! Cbar
         XPQ(-4)= CHM
*! B
         XPQ(5)= BOT
*! Bbar
         XPQ(-5)= BOT
*! T
         XPQ(6)= TOP
*! Tbar
         XPQ(-6)= TOP
      ENDIF

      RETURN
      END
