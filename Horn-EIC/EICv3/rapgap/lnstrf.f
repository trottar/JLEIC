*CMZ :  2.06/39 16/07/98  11.38.10  by  Hannes Jung
*CMZ :  4.41/00 10/10/97  08.39.21  by  Hannes Jung
*-- Author :
C######################################################################
C
C   Various routines to give structure function parametrizations.
C   All but LNSTRF can be used separately without initialization.
C
C ********************************************************************

      SUBROUTINE LNSTRF(X,Q2,XPQ)

C...Structure function per nucleon for a proton/neutron mixture
C...according to defined nucleus.

      DIMENSION XPQ(-6:6)

      CALL PYSTFU(2212,X,Q2,XPQ)

      RETURN
      END
