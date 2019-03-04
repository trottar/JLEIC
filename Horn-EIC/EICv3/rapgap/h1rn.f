*CMZ :  2.08/04 22/12/99  18.07.27  by  Hannes Jung
*CMZ :  2.06/02 09/09/97  11.19.50  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
*-- AUTHOR :  STEPHAN EGLI
      FUNCTION H1RN()
************************************************************************
* Return one random number between 0 and 1 (excl.)                     *
************************************************************************
      DIMENSION X(1)

      CALL H1RNV(X,1)
      H1RN=X(1)

      RETURN
      END
