      SUBROUTINE ANALYS(ntID)

      IMPLICIT NONE

      integer ntID

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      include 'xpikin.inc'

      LOGICAL lEvTyp
      EXTERNAL lEvTyp

      LOGICAL first
      DATA first /.true./

      INTEGER nevent
      DATA nevent /0/
      SAVE nevent

      REAL pi
      PARAMETER (pi = 3.14159265359)

      INTEGER i

      INTEGER iMassOff
      save iMassOff

      REAL hrsAng, bbAng
      PARAMETER (hrsAng = pi-19.6/180.0*pi)
      PARAMETER (bbAng = pi-30.0/180.0*pi)

      real hrsDir(3), bbDir(3)
      save hrsDir, bbDir

      real bbEang, hrsEang, bbNang

      real pHRS, pBBmin, pBBmax
      parameter (pHRS = 2.0)
      parameter (pBBmin = 0.250, 
     &           pBBmax = 0.900)

      logical eInHrs, eInBB, nInBB

C**** Begin ****

      IF (first) THEN
         first = .false.
         
         iMassOff = 0

         hrsDir(1) = sin(hrsAng)
         hrsDir(2) = 0.0
         hrsDir(3) = cos(hrsAng)

         bbDir(1) = sin(bbAng)
         bbDir(2) = 0.0
         bbDir(3) = cos(bbAng)

      ENDIF

      call kinCalc
      if (iNeu .eq. 0) then
C       no neutron/proton found
        return
      endif
      if (abs(amN2 - 0.93956**2) .gt. 0.0025) then
C       problems with these events not matching neutron mass
C       including them in the distribution forms a hyperbolic curve
C       in thetN vs thetE
        iMassOff = iMassOff + 1
        return
      endif

      pDotQ  = - dp(2,1)*dp(iQ,1) - dp(2,2)*dp(iQ,2) - dp(2,3)*dp(iQ,3) 
     &        + dp(2,4)*dp(iQ,4)
      IF (abs(pDotQ) .le. small) pDotQ = sign(small, pDotQ)


      ppDotQ = - dp(iNeu,1)*dp(iQ,1) - dp(iNeu,2)*dp(iQ,2)
     &         - dp(iNeu,3)*dp(iQ,3) + dp(iNeu,4)*dp(iQ,4)

C      IF (.not. lEvTyp()) THEN
C        write(6,100) 
C 100    format(x, 'event rejected')
C        call lulist(1)
C        RETURN
C      ENDIF

      p3Pi_2 = dp(iPi,1)**2 + dp(iPi,2)**2 + dp(iPi,3)**2
      p3Pi = sqrt(p3Pi_2) 

      pTN = sqrt(dp(iNeu,1)**2 + dp(iNeu,2)**2)
      betaN = pN / eN
C     see jackson, classical electrodyanmics
      eNKin = aMassN * (1.0 / sqrt( 1 - betaN**2) - 1.0)

      pTE = sqrt(dp(4,1)**2 + dp(4,2)**2)

C     now find some smearing
      pNOld(1) = dp(iNeu, 1)
      pNOld(2) = dp(iNeu, 2)
      pNOld(3) = dp(iNeu, 3)
      pNOld(4) = dp(iNeu, 4)
      call smearPN(pNNew, pNOld)

      pNNewDotQ =  - pNNew(1)*dp(3,1) - pNNew(2)*dp(3,2) 
     &             - pNNew(3)*dp(3,3) 
     &        + pNNew(4)*dp(3,4)

      xLNew = pNNewDotQ / pDotQ
      IF (abs(1.0 - xLNew) .gt. small) THEN
        xPiNew = x / (1.0 - xLNew)
      ELSE
        xPiNew = x / sign(small, 1.0 - xLNew)
      ENDIF

      call fillHist(100)
      if (xL .ge. 0.75) then
         call fillHist(200)
         call hfnt(ntID)
      endif

      return

      end

      subroutine fillHist(idOff)

      IMPLICIT NONE

      integer idOff

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      include 'xpikin.inc'

C**** Begin ****

      call hf1(idOff + 01, x, 1.0)
      call hf1(idOff + 02, -q2, 1.0)
      call hf2(idOff + 03, x, -q2, 1.0)
      call hf1(idOff + 04, xL, 1.0)
      call hf2(idOff + 05, xL, -q2, 1.0)
      call hf1(idOff + 06, xPi, 1.0)
      call hf2(idOff + 07, xPi, -q2, 1.0)
      call hf2(idOff + 08, xPi, xL, 1.0)
 
      call hf1(idOff + 11, log(x)/ln10, 1.0)
      call hf1(idOff + 12, log(-q2)/ln10, 1.0)
      call hf2(idOff + 13, log(x)/ln10, log(-q2)/ln10, 1.0)
      call hf1(idOff + 14, log(xL)/ln10, 1.0)
      call hf2(idOff + 15, log(xL)/ln10, log(-q2)/ln10, 1.0)
      call hf1(idOff + 16, log(xPi)/ln10, 1.0)
      call hf2(idOff + 17, log(xPi)/ln10, log(-q2)/ln10, 1.0)

      call hf2(idOff + 21, pTN, thetN, 1.0)
      call hf2(idOff + 22, -q2, thetE, 1.0)
      call hf1(idOff + 23, thetN, 1.0)
      call hf1(idOff + 24, thetE, 1.0)
      call hf2(idOff + 25, thetE, thetN, 1.0)
      call hf2(idOff + 26, phiE, phiN, 1.0)
      call hf1(idOff + 27, pTN, 1.0)
      call hf1(idOff + 28, pN, 1.0)
      call hf1(idOff + 29, eN, 1.0)
      call hf2(idOff + 30, pN, thetN, 1.0)
      call hf1(idOff + 31, eNKin, 1.0)

      call hf1(idOff + 41, pdotq, 1.0)
      call hf1(idOff + 42, ppdotq, 1.0)

      call hf2(idOff + 43, amPi2, xPi, 1.0)
      call hf2(idOff + 44, xPi, thetN, 1.0)
      call hf1(idOff + 45, p3Pi, 1.0)
      call hf1(idOff + 46, amPi2, 1.0)
      call hf1(idOff + 47, amE2, 1.0)
      call hf1(idOff + 48, amN2, 1.0)

      call hf2(idOff + 50, xPiNew, xPi, 1.0)
      call hf2(idOff + 51, xLNew, xL, 1.0)

      call hf2(idOff + 55, amPi2, xL, 1.0)
      call hf2(idOff + 56, amPi2, -q2, 1.0)


      RETURN

      END
C
      SUBROUTINE bookHist(idOff)

      IMPLICIT NONE

      INTEGER idOff

C**** Begin ****

      call hbook1(idOff + 01, 'x', 100, 0.0, 1.0, 0.0)
      call hbook1(idOff + 02, 'q^2!', 100, 0.0, 50.0, 0.0)
      call hbook2(idOff + 03, 'q^2! vs x', 50, 0.0, 1.0,
     &                             50, 0.0, 50.0, 0.0)
      call hbook1(idOff + 04, 'x?L!', 100, 0.0, 1.0, 0.0)
      call hbook2(idOff + 05, 'q^2! vs x?L!', 50, 0.0, 1.0,
     &                             50, 0.0, 50.0, 0.0)
      call hbook1(idOff + 06, 'x?[p]!', 100, 0.0, 1.0, 0.0)
      call hbook2(idOff + 07, 'q^2! vs x?[p]!', 50, 0.0, 1.0,
     &                             50, 0.0, 10.0, 0.0)
      call hbook2(idOff + 08, 'x?L! vs x?[p]!', 50, 0.0, 1.0,
     &                             50, 0.0, 1.0, 0.0)

      call hbook1(idOff + 11, 'log(x)', 100, -5.0, 0.0, 0.0)
      call hbook1(idOff + 12, 'log(q^2!)', 100, 0.0, 1.5, 0.0)
      call hbook2(idOff + 13, 'log(q^2!) vx log(x)', 50, -5.0, 0.0,
     &                             50, 0.0, 2.0, 0.0)
      call hbook1(idOff + 14, 'log(x?L!)', 100, -5.0, 0.0, 0.0)
      call hbook2(idOff + 15, 'log(q^2!) vs log(x?L!))', 50, -5.0, 0.0,
     &                             50, 0.0, 1.5, 0.0)
      call hbook1(idOff + 16, 'log(x?[p]!)', 100, -5.0, 0.0, 0.0)
      call hbook2(idOff + 17, 'log(q^2!) vs log(x?[p]!)', 
     &                             50, -5.0, 0.0,
     &                             50, 0.0, 1.0, 0.0)


      call hbook2(idOff + 21, 'p?T! vs [q]?N!',  75, 0.0, 1.0,
     &            75, 0.0, 3.1416, 0.0)
      call hbook2(idOff + 22, '[q]?e! vs Q^2!', 75, 0.0, 50.0,
     &            75, 0.0, 3.1416, 0.0)
      call hbook1(idOff + 23, '[q]?N!', 200, 0.0, 3.1416, 0.0)
      call hbook1(idOff + 24, '[q]?e!', 200, 0.0, 3.1416, 0.0)
      call hbook2(idOff + 25, '[q]?N! vs [q]?e!', 75, 0.0, 3.1416,
     &            75, 0.0, 3.1416, 0.0)
      call hbook2(idOff + 26, '[f]?N! vs [f]?e!', 50, -3.1416, 3.1416,
     &            50, 0.0, 3.1416, 0.0)
      call hbook1(idOff + 27, 'p^N!?T!', 100, 0.0, 1.0, 0.0)
      call hbook1(idOff + 28, 'p?N!', 100, 0.0, 1.0, 0.0)
      call hbook1(idOff + 29, 'e?N!', 100, 0.0, 2.0, 0.0)
      call hbook2(idOff + 30, 'p^N!?T! vs [q]?N!',  75, 0.0, 1.0,
     &            75, 0.0, 3.1416, 0.0)
      call hbook1(idOff + 31, 'E^N!?ken!', 100, 0.0, 1.0, 0.0)

      call hbook1(idOff+41, 'p dot q', 100, 0.0, 50.0, 0.0)
      call hbook1(idOff+42, 'pp dot q', 100, 0.0, 40.0, 0.0)

      call hbook2(idOff+43, 'x?[p]! vs m?[p]!^2', 50, -1.2, 0.0, 
     &            50, 0.0, 1.0, 0.0)
      call hbook2(idOff+44, '[q]?N! vs x?[p]!', 50, 0.0, 1.0,
     &            50, 0.0, 3.1416, 0.0)

      call hbook1(idOff+45, 'abs(p?[p]!)', 100, 0.0, 5.0, 0.0)
      call hbook1(idOff+46, 'm?[p]!^2', 100, -1.2, 0.0, 0.0)
      call hbook1(idOff+47, 'm?[e]!^2', 100, -0.4e-5, 0.4e-5, 0.0)
      call hbook1(idOff+48, 'm?N!^2', 100, 0.8, 1.0, 0.0)

      call hbook2(idOff+50, 'x?[p]! vs x?[p]!^meas!', 50, 0.0, 1.0, 
     &        50, 0.0, 1.0, 0.0)

      call hbook2(idOff+51, 'x?L! vs x?L!^meas!', 50, 0.0, 1.0, 
     &        50, 0.0, 1.0, 0.0)

      call hbook2(idOff+55, 'x?L! vs m?[p]!^2', 50, -1.2, 0.0, 
     &            50, 0.0, 1.0, 0.0)
      call hbook2(idOff+56, 'Q^2! vs m?[p]!^2', 50, -1.2, 0.0, 
     &            50, 0.0, 50.0, 0.0)



      RETURN

      END
C
      SUBROUTINE bookNtuple(ntID)

C     This routine will book an ntuple to hold the contents of the 
C     Lujets common block as well as some other interesting 
C     variables (all of which can be calculated from lujets).

      IMPLICIT none

      integer ntID

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      include 'xpikin.inc'

C**** Begin ****

      call hbnt(ntID, 'rapgap-lujets', 'd')

      call hbname(ntID, 'lujets', n, 'npart')

      call hbname(ntID, 'xvalue', x, 'x')
      call hbname(ntID, 'xvalue', xPi, 'xPi')
      call hbname(ntID, 'xvalue', xPiFree, 'xPiFree')
      call hbname(ntID, 'xvalue', xL, 'xL')

      call hbname(ntID, 'pion', pPi, 'pPi')
      call hbname(ntID, 'pion', ePi, 'ePi')
      call hbname(ntID, 'pion', amPi2, 'amPi2')
      call hbname(ntID, 'pion', pxPi, 'pxPi')
      call hbname(ntID, 'pion', pyPi, 'pyPi')
      call hbname(ntID, 'pion', pzPi, 'pzPi')

      call hbname(ntID, 'photon', q2, 'q2')
      call hbname(ntID, 'photon', pGam, 'pGam')
      call hbname(ntID, 'photon', eGam, 'eGam')
      call hbname(ntID, 'photon', pxGam, 'pxGam')
      call hbname(ntID, 'photon', pyGam, 'pyGam')
      call hbname(ntID, 'photon', pzGam, 'pzGam')

      call hbname(ntID, 'electron', pE, 'pE')
      call hbname(ntID, 'electron', eE, 'eE')
      call hbname(ntID, 'electron', amE2, 'amE2')
      call hbname(ntID, 'electron', pxE, 'pxE')
      call hbname(ntID, 'electron', pyE, 'pyE')
      call hbname(ntID, 'electron', pzE, 'pzE')
      call hbname(ntID, 'electron', thetE, 'thetE')
      call hbname(ntID, 'electron', phiE, 'phiE')

      call hbname(ntID, 'neutron', iNeu, 'iNeu')
      call hbname(ntID, 'neutron', pN, 'pN')
      call hbname(ntID, 'neutron', eN, 'eN')
      call hbname(ntID, 'neutron', amN2, 'amN2')
      call hbname(ntID, 'neutron', pxN, 'pxN')
      call hbname(ntID, 'neutron', pyN, 'pyN')
      call hbname(ntID, 'neutron', pzN, 'pzN')
      call hbname(ntID, 'neutron', thetN, 'thetN')
      call hbname(ntID, 'neutron', phiN, 'phiN')


      RETURN

      END
C
      SUBROUTINE kinCalc

C     This routine will calculate kinematic quanties in common blocks.

      IMPLICIT none

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      REAL amPiReal
      PARAMETER (amPiReal = 0.13957018)

      include 'xpikin.inc'
      
      REAL pPiDotQ
      INTEGER i

      REAL*8 ePi_D, pPi_D, eN_D, pN_D, eE_D, pE_D

C      real q2check

C**** Begin ****

      iNeu = 0
      DO i = 5, n
        IF (((k(i,2) .eq. 2212) .and. (k(i,4) .eq. 0)) .or.
     &      ((k(i,2) .eq. 2112) .and. (k(i,4) .eq. 0))) iNeu = i
      ENDDO

      IF (iNeu .eq. 0) THEN
C       no neutron/proton found
C        write(6,*) 'No neutron/proton found'
        return
      endif

      iE = 4
      iPi = 5
      iQ = 3

      pPi_D = dsqrt(dp(iPi,1)**2 + dp(iPi,2)**2 + dp(iPi,3)**2)
      pPi = sngl(pPi_D)
      ePi_D = dp(iPi,4)
      ePi = p(iPi,4)
      amPi2 = sngl((ePi_D - pPi_D) * (ePi_D + pPi_D))
      pxPi = p(iPi,1)
      pyPi = p(iPi,2)
      pzPi = p(iPi,3)

      pGam = dsqrt(dp(iQ,1)**2 + dp(iQ,2)**2 + dp(iQ,3)**2)
      eGam = p(iQ,4)
      pxGam = p(iQ,1)
      pyGam = p(iQ,2)
      pzGam = p(iQ,3)

      thetE = datan2(dsqrt(dp(iE,1)**2 + dp(iE,2)**2), dp(iE,3))
      phiE = datan2(dp(iE,2), dp(iE,1))
      pE_D = dsqrt(dp(iE,1)**2 + dp(iE,2)**2 + dp(iE,3)**2)
      pE = sngl(pE_D)
      eE_D = dp(iE,4)
      eE = p(iE,4)
      amE2 = sngl((eE_D - pE_D) * (eE_D + pE_D))
      pxE = p(iE,1)
      pyE = p(iE,2)
      pzE = p(iE,3)
      eDir(1) = p(iE,1) / pE_D
      eDir(2) = p(iE,2) / pE_D
      eDir(3) = p(iE,3) / pE_D

      thetN = datan2(dsqrt(dp(iNeu,1)**2 + dp(iNeu,2)**2), dp(iNeu,3))
      phiN = datan2(dp(iNeu,2), dp(iNeu,1))
      pN_D = dsqrt(dp(iNeu,1)**2 + dp(iNeu,2)**2 + dp(iNeu,3)**2)
      pN = sngl(pN_D)
      eN_D = dp(iNeu,4)
      eN = p(iNeu,4)
      amN2 = sngl((eN_D - pN_D) * (eN_D + pN_D))
      pxN = p(iNeu,1)
      pyN = p(iNeu,2)
      pzN = p(iNeu,3)
      nDir(1) = p(iNeu,1) / pN_D
      nDir(2) = p(iNeu,2) / pN_D
      nDir(3) = p(iNeu,3) / pN_D

      pDotQ  = - dp(2,1)*dp(iQ,1) - dp(2,2)*dp(iQ,2) 
     &         - dp(2,3)*dp(iQ,3) + dp(2,4)*dp(iQ,4)
      IF (abs(pDotQ) .le. small) pDotQ = sign(small, pDotQ)

      pPiDotQ  = - dp(iPi,1)*dp(iQ,1) - dp(iPi,2)*dp(iQ,2)
     &           - dp(iPi,3)*dp(iQ,3) + dp(iPi,4)*dp(iQ,4)
      IF (abs(pPiDotQ) .le. small) pPiDotQ = sign(small, pPiDotQ)

      q2 =  - dp(iQ,1)*dp(iQ,1) - dp(iQ,2)*dp(iQ,2)
     &      - dp(iQ,3)*dp(iQ,3) + dp(iQ,4)*dp(iQ,4)
C      q2check = (egam - pgam) * (egam + pgam)
C      if (abs(q2check - q2) .ge. small) then
C         write(6,*) 'q2check', q2, q2check
C      endif
      x =  -q2 / (2.0 * pDotQ)

      ppDotQ = - dp(iNeu,1)*dp(iQ,1) - dp(iNeu,2)*dp(iQ,2)
     &         - dp(iNeu,3)*dp(iQ,3) + dp(iNeu,4)*dp(iQ,4)
      xL = ppDotQ / pDotQ
      IF (abs(1.0 - xL) .gt. small) THEN
        xPi = x / (1.0 - xL)
      ELSE
        xPi = x / sign(small, 1.0 - xL)
      ENDIF

      xPiFree = -q2 / (2.0 * pPiDotQ)

      RETURN

      END
C
      LOGICAL FUNCTION lEvTyp()

C     Check to see if this is the correct event type

      IMPLICIT NONE

      INTEGER LUPAN
      PARAMETER (LUPAN=4000)

      INTEGER N,K
      REAL P,V
      DOUBLE PRECISION DP
      COMMON/LUJETS/N,K(LUPAN,5),P(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/DP(LUPAN,5)

      REAL tolerance
      PARAMETER (tolerance = 1.0e-3)

      INTEGER i

C**** Begin ****

      lEvTyp = .false.

      IF (n .lt. 4) THEN
C       need 4 particles in record (e, e', gamma, N)
        call reject(1)
        RETURN
      ENDIF
      IF ((abs(k(1,2)) .ne. 11) .or. (abs(k(2,2)) .ne. 2212) .or.
     &    (abs(k(4,2)) .ne. 11) ) THEN
C       first three particles need to be initial electron, proton 
C       and scattered electron
        call reject(2)
        RETURN
      ENDIF
      IF ((abs(dp(4,1)+dp(3,1)-dp(1,1)) .ge. tolerance) .or.
     &    (abs(dp(4,2)+dp(3,2)-dp(1,2)) .ge. tolerance) .or.
     &    (abs(dp(4,3)+dp(3,3)-dp(1,3)) .ge. tolerance) .or.
     &    (abs(dp(4,4)+dp(3,4)-dp(1,4)) .ge. tolerance) ) THEN
C       Check for momentum balance in e e' vertex--ie this is an 
C       elastic scattering event
C        write(6,*) 1, dp(4,1)+dp(3,1)-dp(1,1)
C        write(6,*) 2, dp(4,2)+dp(3,2)-dp(1,2)
C        write(6,*) 3, dp(4,3)+dp(3,3)-dp(1,3)
C        write(6,*) 4, dp(4,4)+dp(3,4)-dp(1,4)
        call reject(3)
        RETURN
      ENDIF


      lEvTyp = .true.

      RETURN

      END
C     
      SUBROUTINE initHbook(nhbmem, ntID)

      implicit none
      
      integer nhbmem
      integer ntID

      CHARACTER ntFileName*80, fileName*70
      INTEGER lenocc
      EXTERNAL lenocc

      Integer istat
      integer irecl
      parameter (irecl = 1024)

      INTEGER iquest(100)
      COMMON /quest/ iquest

C**** Begin ****

      filename = ' '
      call getenv('RGFIXED_EIC', filename)
      if (filename(1:1) .eq. ' ') then
        filename = 'rapgap'
      endif

      call hlimit(nhbmem)

      ntFileName = filename(1:lenocc(filename))//'.nt'
      iquest(10) =  64000
      call hropen(10, 'rapgapNT', ntFileName, 'n', irecl, istat)
      call bookNtuple(ntID)

      call hcdir('//PAWC', ' ')
      call bookHist(100)
      call bookHist(200)

      return

      end
C
      Subroutine closeHbook(ntID)

      implicit none
      
      integer ntID

      CHARACTER histFileName*80, fileName*70
      INTEGER lenocc
      EXTERNAL lenocc

      integer icycle, istat
      integer irecl
      parameter (irecl = 1024)

C**** Begin ****
      

      filename = ' '
      call getenv('RGFIXED_EIC', filename)
      if (filename(1:1) .eq. ' ') then
        filename = 'rapgap'
      endif
      
C      close ntuple
      call hcdir('//rapgapNT', ' ')
      call hrout(ntID, icycle, ' ')
      call hrendc('rapgapNT')
      close(10)
      call hcdir('//PAWC', ' ')
      call hdelet(ntID)

      histFileName = filename(1:lenocc(filename))//'.his'
      call hropen(10, 'rapgap', histFileName, 'n', irecl, istat)
      call hrout(0, icycle, ' ')
      call hrendc('rapgap')
      close(10)
      
      return

      end
C
      
      SUBROUTINE reject(iReason)

      IMPLICIT NONE

      INTEGER iReason

      INTEGER i

      INTEGER mxReason
      PARAMETER (mxReason = 10)
      INTEGER nReject(mxReason)
      SAVE nReject

      LOGICAL first
      DATA first /.true./
      SAVE first

C**** Begin ****

      IF (first) THEN
        DO i = 1, mxReason
          nReject(i) = 0
        ENDDO
        first = .false.
      ENDIF

      IF ((0 .lt. iReason) .and. (iReason .le. mxReason)) THEN
        nReject(iReason) = nReject(iReason) + 1
      ELSE IF (iReason .eq. -1) THEN
C       signal to print rejections
        write(6,100)
 100    format(x, '*********************************', /, 
     &         'Summary of reasons for rejecting event:')
        DO i = 1, mxReason
          write(6,101) i, nReject(i)
 101      format(x, i3, i6)
        ENDDO
        write(6,102) 
 102    format(x, 'End of rejection summary', /
     &         '*********************************')
      ELSE
       write(6,103), iReason
 103   format(x, 'Rejection reason', i3, ' not defined.')
      ENDIF

      RETURN

      END
C
      SUBROUTINE smearPN(pNNew, p)

      IMPLICIT none

      REAL pNNew(4), p(4)
      REAL pTot, pTotOld

      REAL pi
      PARAMETER (pi = 3.14159265358979323846264338327950288)

      REAL sigAngle, sigEnergy
      PARAMETER (sigAngle = 2.5 / 180.0 * pi) ! n degree resolution
      PARAMETER (sigEnergy = 0.06) ! 5% measurement

      real dEnergy, alpha
      REAL aMass2, energy

      REAL sinAlpha, cosAlpha, gamma
      REAL cosph, costh, sinph, sinth
      REAL v(3), s(3), sNorm

C****************************
C**** note that since this is linking with the librapgap.a prior to the
C**** link to cernlib, rndm is taken from librapgap.a, and this is a 
C**** Double precision function
C****************************
C      REAL rndm, rdummy
C      EXTERNAL rndm

      REAL h1rn
      EXTERNAL h1rn

C**** Begin ****

      call rannor(alpha, dEnergy)
      alpha = alpha * sigAngle
      sinAlpha = dble(sin(alpha))
      cosAlpha = dble(cos(alpha))

C      gamma = 2.0 * pi * rndm(rdummy)
      gamma = 2.0 * pi * h1rn()

C     new direction vector, relative to old direction
C     aligned along v(3) axis
      v(1) = sinAlpha * cos(gamma)
      v(2) = sinAlpha * sin(gamma)
      v(3) = cosAlpha

C     Calculate rotation matrix back to lab axis.  First
C     rotate around z axis by ph:
C
C             / cos(ph) -sin(ph)  0 \
C     R(ph) = | sin(ph)  cos(ph)  0 |
C             \   0        0      1 /
C
C     next rotate around y axiz by th.
C
C             /  cos(th) 0 sin(th) \
C     R(th) = |    0     1   0     |
C             \ -sin(th) 0 cos(th) /
C
C     and the product of these two rotations:
C
C     R(ph) * R(th) = 
C
C         /  cos(th)*cos(ph) -sin(ph) sin(th)*cos(ph) \
C         |  cos(th)*sin(ph)  cos(ph) sin(th)*sin(ph) |
C         \     -sin(th)        0         cos(th)     /
C

      pTot = sqrt(p(1)**2+p(2)**2+p(3)**2)

      IF ((abs(p(1)) .gt. 1e-10) .or.
     &    (abs(p(2)) .gt. 1e-10)) THEN
        costh = p(3) / pTot
        sinth = sqrt(p(1)**2+p(2)**2) / pTot
        cosph = p(1) / sqrt(p(1)**2+p(2)**2)
        sinph = p(2) / sqrt(p(1)**2+p(2)**2)
      ELSE
        costh = 1.0
        sinth = 0.0
        cosph = 1.0
        sinph = 0.0
      ENDIF

      s(1) = costh*cosph*v(1) - sinph*v(2) + sinth*cosph*v(3)
      s(2) = costh*sinph*v(1) + cosph*v(2) + sinth*sinph*v(3)
      s(3) =   -sinth*v(1)    +                 costh*v(3)

      sNorm = sqrt(s(1)**2 + s(2)**2 + s(3)**2)

      aMass2 = p(4)**2 - p(1)**2 - p(2)**2 - p(3)**2
      pTotOld = pTot
      pTot = max(0.0d0, pTot * dble(1.0 + dEnergy * sigEnergy) )
      energy = sqrt(aMass2 + pTot**2)

C      pNNew(1) = p(1)
C      pNNew(2) = p(2)
C      pNNew(3) = p(3)
C      pNNew(4) = p(4)
      pNNew(1) = s(1) / sNorm * pTot
      pNNew(2) = s(2) / sNorm * pTot
      pNNew(3) = s(3) / sNorm * pTot
      pNNew(4) = energy

      RETURN

      END
       
