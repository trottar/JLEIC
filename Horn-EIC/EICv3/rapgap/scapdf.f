*CMZ :  2.08/05 28/03/2000  14.07.14  by  Hannes Jung
*-- Author :    Hannes Jung   25/04/96
      SUBROUTINE SCAPDF(BETA,SCALE,XPQ,X_POM,T2)

      Implicit None

      REAL T2,X_POM
      REAL BETA,SCALE
      REAL XPQ(-6:6)
      integer i

      DO I=-6,6
         XPQ(I)=0.0
      ENDDO

      CALL SCA_MAIN(beta,scale,xpq,x_pom,t2)

      RETURN

      END
