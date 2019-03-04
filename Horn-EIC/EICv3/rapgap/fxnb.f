*CMZ :  2.08/00 06/06/99  17.40.09  by  Hannes Jung
*CMZ :  2.07/03 24/05/99  17.32.18  by  Hannes Jung
*-- Author :    Hannes Jung   03/04/94
      FUNCTION FXNB(X)
      IMPLICIT None
	Double Precision FXNB,FXNB1,X,XG,XGF,FXN1,WEIGHT
	Integer MXDIM,I
      PARAMETER (MXDIM = 50)
      DIMENSION X(MXDIM),XGF(20)
      COMMON/XFXNB/XG(20)
      EXTERNAL FXN1
      FXNB = 0.0D0
	Weight = 0.D0
c      write(6,*) ' fxnb ',(x(i),i=1,10)
      DO 10 I=1,20
         XGF(I) = X(I)
         XG(I) = X(I)
   10 CONTINUE
c      write(6,*) ' fxnb ',xgf
      FXNB1 = FXN1(XGF,WEIGHT)
c      write(6,*) ' fxnb ',fxnb1
      FXNB = FXNB1
      RETURN
      END
