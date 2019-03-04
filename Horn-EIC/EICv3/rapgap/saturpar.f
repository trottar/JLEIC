*CMZ :  2.08/05 27/03/2000  15.56.08  by  Hannes Jung
*-- Author :    Hannes Jung   27/03/2000
      Subroutine SATURPAR
      Implicit NONE

      Double Precision  R02, LAM, SIGM0, X0, ALPHAS
      COMMON /SATURM/   R02, LAM, SIGM0, X0, ALPHAS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Parameters of Wuesthoff Model
C Alphas
           ALPHAS = 0.2
C SIGM0 is 23.0 mB we write this in GeV-2 as 23.03*2.568
c ori      SIGM0 = 23.03*2.568
           SIGM0 = 23.03*2.568
c ori      LAM  = 0.288
           LAM  = 0.288
c ori      X0 =  3.04/10000.
           X0 =  3.04/10000.

      RETURN
      END
