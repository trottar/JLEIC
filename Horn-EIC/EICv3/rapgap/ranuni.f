
*CMZ :  2.08/05 27/03/2000  16.15.50  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      FUNCTION RANUNI(A,B)
C ... Uniform random random number in range [A,B]
      DOUBLE PRECISION RANUNI,A,B,RN
      Real RVec(5)
	Double Precision draprn

ccc      Call RANMAR(RVec,1)
c      RN=RVec(1)
c changed to standard draprn

      RN=draprn()
      RANUNI=A+RN*(B-A)
      END
