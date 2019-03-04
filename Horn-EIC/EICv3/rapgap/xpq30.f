*CMZ :  2.08/05 27/03/2000  16.05.12  by  Hannes Jung
*CMZ :  2.08/04 09/01/2000  13.21.10  by  Hannes Jung
*CMZ :  2.08/02 14/07/99  18.39.56  by  Hannes Jung
*-- Author :    Hannes Jung   14/07/99
      Subroutine XPQ30(X,Q2,XPQ)
      Implicit None
      Integer I
      REAL X,Q2,XPQ
      Double Precision XD,Q2D,SIGOUT,ESIGOUT,WEIOUT
      Double Precision pdfnorm
      Dimension XPQ(-6:6)
*KEEP,RGPARAM.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
*KEND.
      DOUBLE PRECISION POM,PIM
      COMMON/WEIGHT1/ POM,PIM
      Double Precision GEV2NB
      DATA GEV2NB/.3893857D+3/

      XD = X
      Q2D = Q2
      Do I=-6,6
         XPQ(I) = 0.
      Enddo

      Call SATSIGTOT(Q2D,XD,SIGOUT,ESIGOUT,WEIOUT)
c     write(6,*) 'XPQ30:Q2D,XD,SIGOUT,ESIGOUT,WEIOUT',
c     &  Q2D,XD,SIGOUT,ESIGOUT,WEIOUT
c factor 9 comes from 1/eq**2 for d_quark with eq=1/3
ccc   XPQ(1) = SIGOUT*Q2D*9.D0/4.D0/PI/PI/ALPH
      pdfnorm=1./2.d0/9.d0
      pdfnorm=Q2D*9.D0/12.D0/4.D0/PI/PI/ALPH/GEV2NB
      XPQ(1) = SIGOUT*pdfnorm
      XPQ(-1) = SIGOUT*pdfnorm
      XPQ(2) = SIGOUT*pdfnorm
      XPQ(-2) = SIGOUT*pdfnorm
      XPQ(3) = SIGOUT*pdfnorm
      XPQ(-3) = SIGOUT*pdfnorm
      POM = WEIOUT
c      write(6,*) 'XPQ30: XPQ(1) ',XPQ(1),' weight ',weiout,' x,q2 ',x,q2

      RETURN
      END
