*CMZ :  2.08/04 22/12/99  15.42.05  by  Hannes Jung
*CMZ :  2.07/03 15/05/99  21.50.10  by  Hannes Jung
*-- Author :
      SUBROUTINE LPRIKT(S,PT,PHI)
	IMPLICIT NONE
      Double Precision draprn
	Real S,PT,PHI
      EXTERNAL draprn
C...Size (PT) and azimuthal angle (PHI) of primordial kt according
C...to a Gaussian distribution.

      PT=S*SQRT(-ALOG(Real(draprn())))
      PHI=6.2832*draprn()
      RETURN
      END
