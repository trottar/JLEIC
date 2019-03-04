*CMZ :  2.08/04 22/12/99  15.39.28  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  17.33.19  by  Hannes Jung
*-- Author :    Hannes Jung   03/05/98
c get draprn for bases/spring
      FUNCTION DRN(ISEED)
	Implicit None
	Integer ISEED
      double precision drn1,DRN
      Double Precision draprn
	EXTERNAL draprn
      DRN1= draprn()
      DRN = DRN1
	RETURN
	END
