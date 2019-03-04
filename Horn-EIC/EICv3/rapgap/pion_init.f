*CMZ :  2.06/29 15/03/98  13.45.09  by  Hannes Jung
*CMZ :  2.06/17 12/08/97  20.38.58  by  Julian Phillips
*CMZ :  2.03/12 11/10/96  13.36.58  by  Hannes Jung
*CMZ :  2.03/08 13/09/96  20.13.30  by  Julian Phillips
*CMZ :  1.00/03 08/07/96  16.48.11  by  Julian Phillips
*-- Author :    Julian Phillips   08/07/96
      SUBROUTINE PION_INIT(NPTYPE,NGROUP,NSET)

C Initialise PDFlib for extraction of pion and creation of meson
      IMPLICIT REAL*8 (A-G,O-Z)
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20),NPTYPE,NGROUP,NSET
      LOGICAL FIRST
      LOGICAL PDFFIRST
      COMMON/W50516/PDFFIRST

      DATA FIRST/.TRUE./
      PDFFIRST = FIRST
C      INTEGER NPTYPE,NGROUP,NSET

C Initialise for Pion Structure Function
      PARM(1) = 'NPTYPE'
      VAL(1) = NPTYPE
      PARM(2) = 'NGROUP'
      VAL(2) = NGROUP
      PARM(3) = 'NSET'
      VAL(3) = NSET

C Inform User of Action
      IF(FIRST) THEN
         FIRST = .FALSE.
         WRITE(6,*)'Initilising for PION structure function'
         WRITE(6,*)PARM(1),VAL(1)
         WRITE(6,*)PARM(2),VAL(2)
         WRITE(6,*)PARM(3),VAL(3)
      ENDIF

      CALL PDFSET(PARM,VAL)

      RETURN
      END
