*CMZ :  2.08/04 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10.53.28  *-- Author :
C================================================================
C================================================================

      Subroutine EVLRD (NU, HEADER, IRR)
C
C     Before calling this routine, the calling program must open an external
C              file for data transfer, provide a filename for the file, and
C              establish the equivalence of that file with unit=Nu, such as:
C
C       Open (Nu, File='FILENAME', Form='UNFORMATTED', status='OLD', Err=... )
C
C       where the file named 'FILENAME' must exist and it has to be written
C       originally by EVLWT.  (cf. above)
C
C       Input parameter:   Nu = unit number for the external file to be
C                               written to. (established in calling program)
C
C                      Header = informational header for the file
C
C       Output parameter: Irr
C
C                       ----------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      CHARACTER HEADER*78, Line*80
C
      PARAMETER (MXX = 1000, MXQ = 100, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX)
      COMMON / PEVLD1 / KF, Nelmt

      Read  (Nu, '(A)') Header

      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, Am1, Am2, Am3, AM4, AM5, AM6

      Read  (Nu, '(A)') Line
      Read  (Nu, *) IPD0, IHDN, IKNL, NfMx, KF, Nelmt

      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, JT,  NG, NTL(NG+1)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) (NTL(I), NTN(I), TLN(I), DTN(I), I =1, NG)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QV(I), TV(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, XCR,  (XV(I), I =1, NX)
C
C                  Since quark = anti-quark for nfl>2 at this stage,
C                  we Read  out only the non-redundent data points
C                  No of flavors = NfMx sea + 1 gluon + 2 valence
      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      CALL RdUpd (UPD, Npts, NU, IRR)

      CLOSE (NU)
      IF (IRR .NE. 0) PRINT *, 'Read error in EVLRD'

      If (NfMx .GE. 3) Then
C                                       Refill Upd for s -> t
         Do 20 Nflv = 3, NfMx
            J0 = (-Nflv + NfMx) * Nblk
C offset = NfMx sea + 1 gluon + (Nflv-1) less-massive quarks = NfMx+Nflv
            J1 = ( Nflv + NfMx) * Nblk
            Do 10 I = 1, Nblk
               Upd (J1 + I) = Upd (J0 + I)
   10       Continue
   20    Continue
      EndIf

      End
