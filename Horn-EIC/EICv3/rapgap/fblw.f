
      FUNCTION FBLW(X)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DOUBLE PRECISION PT2GEN,PHIGEN
      COMMON/HARDPOM/PT2GEN,PHIGEN
      common   /phi/       phi
      common   /parameter/ s,q,beta,xp,t
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGPARAS.
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE

*KEND.
      REAL * 4 X,FBLW,ULMASS
      EXTERNAL ULMASS
      aem=1.d0/137.d0
      pi=2.d0*dasin(1.d0)
      pt2min=pt2cut(12)
      xm=q/beta-q + t
      pt2max=xm/4.d0 - 4.d0*DBLE(ULMASS(3))**2
      F2QT = 0.D0

      FBLW = 0.0
c      write(6,*) ' FBLW ',pt2min,pt2max,x
      if(pt2min.gt.pt2max) return
      rn1 = dble(x)
c weighting with 1/pt**2
c      pt2 = pt2min*(pt2max/pt2min)**rn1
c      wpt2 = log(pt2max/pt2min)*pt2
c weighting with 1/pt**4
      pt2 = pt2max*pt2min/(pt2max + rn1*(pt2min-pt2max))
      wpt2 = pt2**2*(pt2max-pt2min)/pt2max/pt2min
      CALL SIGBLW(pt2,sigt,sigl)
c         write(6,*) 'f2blw sig,pt2 ',f2qt,pt2
      F2QT=sigt*q/(4.D0*PI**2*aem)*wpt2
      FLQ=sigL*q/(4.D0*PI**2*aem)*wpt2
c          write(6,*) ' f2blw xng,xpom',,xng,xpom
c       write(6,*) ' s,q,beta,xp,t ',s,q,beta,xp,t
c       write(6,*) 'f2blw f2qt ',f2qt,pt2,x
      FBLW = SNGL(F2QT)
      RETURN
      END
