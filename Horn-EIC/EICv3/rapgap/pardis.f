*CMZ :  2.08/04 14/07/99  18.40.56  by  Hannes Jung
*CMZ :  2.07/02 23/03/99  16.46.55  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10*-- Author :
      FUNCTION PARDIS (IPRTN, XX, QQ)
C
C       Given the parton distribution function in the array U in
C       COMMON / PEVLDT / , this routine fetches u(fl, x, q) at any value of
C       x and q using Mth-order polinomial (II=0) or rational fraction (II=1)
C       interpolation for x. It always uses quadratic polinomial interpolation
C       in ln ln (Q/lambda).

C       The calling program must ensure that 0 =< x =< 1 ;
C       If 0 =< x < Xmin, extrapolation is used and a warning is given

C       The calling program must ensure that Alambda < Q ;
C       If (Alambda < Q < Qini .or. Qmax < Q),
C          extrapolation is used and a warning is given
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Character Msg*80
      LOGICAL LSTX
C
      PARAMETER (MXX = 1000, MXQ = 100, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (Smll = 1D-9)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX)
      COMMON / PEVLD1 / KF, Nelmt
C
      Dimension Fq(5), Df(5)
C     Dimension QLg(5)
C                           M determines the order of the polynomial
C            II switches between the polint/ratint interpolation routines
C                                     They are fixed for this version
      Save
      Data M , II / 2, 0 /
      Data Iwrn1, Iwarn2, Iwarn3 / 3*0 /
C
C                                             Check integrity of Upd
C      If (Irun .Eq. 0) Then
C         Print '(1pE13.5, 5E13.5)', (UPD(I), I=1,Nelmt-1)
C         Irun = 1
C      Endif

      X = XX
      Q = QQ
C      M = Max (M, 1)
C      M = Min (M, 3)
      Md = M / 2
      Amd= Md
C
      IF (X .LT. Xmin-Smll) THEN
         Msg = '0 < X < Xmin in ParDis; extrapolation used!'
         CALL WARNR (IWRN3, NWRT, Msg, 'X', X, Xmin, 1D0, 1)
      EndIf
C
      IF (Q .LT. QINI-Smll) THEN
         Msg = 'Q less than QINI in PARDIS call; Q SET = QINI.'
         CALL WARNR (IWRN1, NWRT, Msg, 'Q', Q, QINI, QMAX, 1)
         Q = QINI
      ElseIF (Q .GT. QMAX) THEN
         Msg = 'Q greater than QMAX in PARDIS call; '
         Msg = Msg // 'Extrapolation will be used.'
         CALL WARNR(IWRN2, NWRT, Msg, 'Q', Q, QINI, QMAX, 1)
      EndIf
C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
   10 If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 10
      Endif

      Jx = JL - (M-1)/2
      If     (Jx .LT. 0) Then
         Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
   20 If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (Q .GT. QV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 20
      Endif

      Jq = JL - (M-1)/2
      If     (Jq .LT. 0) Then
         Jq = 0
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
      Endif

      SS  = LOG (Q/AL)
      TT  = LOG (SS)
C                             Find the off-set in the linear array Upd
      JFL = IPRTN + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                    Nt  -| ..........................
C                         | ..........................
C                   Jq+M -| .....o ......o............  Iq=M+1
C                         | .........X................
C                    Jq  -| .(J0)o ......o ...........  Iq=1
C                         | ...... ...................
C                     0  --------|-------|-----------|
C                         0     Jx     Jx+M          Nx

C      Write (Nwrt, '(/ 10(1pE11.3)/)') X, (XV(I),I=Jx,Jx+M), Q

      Do 30 Iq = 1, M+1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
C                                                     Interpolate in x
         If (II .EQ. 0) Then
            Call Polint (XV(Jx), Upd(J1), M+1, X, Fq(Iq), Df(Iq))
         Elseif (II .EQ. 1) Then
            Call RatInt (XV(Jx), Upd(J1), M+1, X, Fq(Iq), Df(Iq))
         Else
            Print *, 'II out of range in Pardis; II = ', II
         Endif

C     QLg(Iq) = Log ( QV(Jq+Iq-1)/AL )
C     Write (Nwrt, '(10(1pE11.3))')
C    >QV(Jq+Iq-1), (Upd(I), I=J1,J1+M), Fq(Iq), Df(Iq)
   30 Continue
C                                                     Interpolate in LnLnQ
      Call Polint (TV(Jq), Fq(1), M+1, TT, Ftmp, Ddf)
C                                                     Or interpolate in LnQ
C     Call Polint (QLg(1), Fq(1), M+1, SS, Ftmp, Ddf)
C     Write (Nwrt, '(/ 10(1pE11.3))') (Fq(I), I=1,M+1), Ftmp, Ddf
C     Write (Nwrt, *) '----------'
C
      PARDIS = Ftmp
C
      RETURN
C                        ****************************
      END
