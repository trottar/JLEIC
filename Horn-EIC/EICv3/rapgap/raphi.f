*CMZ :  2.07/03 24/05/99  17.34.39  by  Hannes Jung
*-- Author :    Hannes Jung   08/05/98
      function raphi(i1,i2,i3)
c calculate azimuthal phi angle
c i1 = electron momentum
c i2 = incoming parton momentum
c i3 = final quark momentum
c cos_phi = (i1 x i2) * (i2 x i3)/abs(i1 x i2)/abs(i2 x i3)
      Implicit none
      Integer i1,i2,i3
      Double precision te11,te12,te13,te21,te22,te23,prod,phi,raphi
*KEEP,RGLUJETS.
      INTEGER N,K
      REAL SP,V
      DOUBLE PRECISION P
      INTEGER LUPAN
      PARAMETER (LUPAN=4000)


      COMMON/LUJETS/N,K(LUPAN,5),SP(LUPAN,5),V(LUPAN,5)
      COMMON/DUJETS/P(LUPAN,5)
      REAL ULMASS
      DOUBLE PRECISION DOT1,DPLU,DLANGL
      EXTERNAL DLANGL,ULMASS,DPLU,LUCHGE,DOT1
C      SAVE

*KEND.
      te11 = P(i2,2)*p(i1,3)-p(i2,3)*p(i1,2)
      te12 = P(i2,3)*p(i1,1)-p(i2,1)*p(i1,3)
      te13 = P(i2,1)*p(i1,2)-p(i2,2)*p(i1,1)
      te21 = P(i2,2)*p(i3,3)-p(i2,3)*p(i3,2)
      te22 = P(i2,3)*p(i3,1)-p(i2,1)*p(i3,3)
      te23 = P(i2,1)*p(i3,2)-p(i2,2)*p(i3,1)
      prod = te11*te21 + te12*te22 + te13*te23
      prod=prod/dsqrt(te11**2+te12**2+te13**2)
      prod=prod/dsqrt(te21**2+te22**2+te23**2)
      if(abs(prod).gt.1.0d0) then
         prod= 1.d0
      endif
      phi = dacos(prod)
      if(phi.ne.phi) then
         write(6,*) ' raphi ',phi,prod
      endif
      raphi = phi
      return
      end
