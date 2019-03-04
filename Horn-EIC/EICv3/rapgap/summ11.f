      FUNCTION SUMM11(X,PHI)
      implicit none
*KEEP,RGPARAM.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
*KEND.
      REAL SUMM11,X,PHI
      DOUBLE PRECISION TIL,F_G,Y,WT
      DOUBLE PRECISION PK1,PK2,PK3,PL
      COMMON /QQG/ PK1(4),PK2(4),PK3(4),PL(4)
      REAL L2MIN,L2MAX
      COMMON /QQG_C/L2MIN,L2MAX
      EXTERNAL TIL,F_G
      Y = DBLE(L2MIN*(L2MAX/L2MIN)**X)
      WT = Y *DBLE(LOG(L2MAX/L2MIN))
      WT = WT/PI
c  phi dep
      PL(1) = DSQRT(Y)*DBLE(COS(PHI))
      PL(2) = DSQRT(Y)*DBLE(SIN(PHI))
      SUMM11 = sngl(WT*TIL(1,1)*F_G(Y)/Y)
c          write(6,*) ' summ11 x,phi ',summ11,x,phi
      IF(SUMM11.LT.0.) THEN
c         write(6,*) ' summ11 x,phi ',x,phi
c         write(6,*) 'TIL(1,1),F_G(Y) ',TIL(1,1),F_G(Y)
      ENDIF
      RETURN
      END
