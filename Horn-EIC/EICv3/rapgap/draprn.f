*CMZ :  2.08/04 22/12/99  18.08.49  by  Hannes Jung
*CMZ :  2.06/02 09/09/97  11.19.50  by  Hannes Jung
*CMZ :  2.04/00 23/12/96  11.44.28  by  Hannes Jung
*CMZ :  1.03/01 03/04/94  16.45.55  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
************************************************************************
      DOUBLE PRECISION FUNCTION draprn ()
************************************************************************
*
*  Purpose: random number with DOUBLE PRECISION in the interval (0,1)
*  -------  based on H1 single precision random generator H1RNV
*
*  Input  : no
*  -----
*
*  Output:  DRAPRN  -- a random number (DOUBLE PRECISION)
*  ------
*
*           Method: two random words with 24-bit mantissa are combined
*                   into one word of 48-bit mantissa
*           Result: dpraprn * 2*48 is a random integer
*                   in the interval (1,2**48-1)
*
*  Author:  A. Fedotov 16.04.1999
*
************************************************************************

      DOUBLE PRECISION TWOM24
      PARAMETER (TWOM24 = 2.D0**(-24))

      DIMENSION R(2)

*-------------------------

*  get two random numbers with single precision
      CALL H1RNV (R, 2)

C      PRINT 900, R(1) * 2.**24, R(2) * 2.**24
C 900  FORMAT (' H1RNDP:', 2F20.5)

*  combine them into one:
      DRAPRN = DBLE (R(1)) + DBLE (R(2)) * TWOM24

      RETURN
      END
