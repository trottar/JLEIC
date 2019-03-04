

      REAL x, xPi, xPiFree, xL
      COMMON /xvalue/ x, xPi, xPiFree, xL

      INTEGER iPi
      REAL ePi, pPi, amPi2, pxPi, pyPi, pzPi
      REAL p3Pi_2, p3Pi
      COMMON /pion/ iPi, ePi, pPi, amPi2, pxPi, pyPi, pzPi, p3pi_2, p3Pi

      INTEGER iQ
      REAL eGam, pGam, pxGam, pyGam, pzGam, q2
      COMMON /photon/ iQ, eGam, pGam, pxGam, pyGam, pzGam, q2

      INTEGER iE
      REAL thetE, phiE, eE, pE, amE2, pxE, pyE, pzE, eDir(3)
      COMMON /electron/ iE, thetE, phiE, eE, pE, amE2, pxE, pyE, pzE,
     &       eDir

      INTEGER iNeu
      REAL thetN, phiN, eN, pN, amN2, pxN, pyN, pzN, nDir(3)

      REAL pTN, pTE
      REAL betaN, eNKin, aMassN
      PARAMETER (aMassN = 0.931494013)

      REAL pNNew(4), pNOld(4)
      REAL xPiNew, xLNew, pNNewDotQ

      COMMON /neutron/ iNeu, thetN, phiN, eN, pN, amN2, pxN, pyN, pzN,
     &       nDir, pTN, pTE, betaN, eNKin, pNNew, pNOld, xPiNew, xLNew,
     &       pNNewDotQ


      REAL ppDotQ, pDotQ
      common /forDots/ ppDotQ, pDotQ

      REAL small
      PARAMETER (small = 1.0e-10)

      REAL ln10
      PARAMETER (ln10 = 2.302585092994)
