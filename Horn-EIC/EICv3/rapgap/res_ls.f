      SUBROUTINE RES_LS(X,Q2,XPQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL XPQ(-6:6),X,Q2
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.
      DIMENSION XX(0:100),Q2X(0:100),XPD(0:100,0:100,-6:6)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      DATA NPR/0/
      IF(FIRST) THEN
         IF(NG.EQ.15) THEN
*         open(30,FILE='diff_tabnew.param.txt', FORM='formatted',STATUS=
            open(30,FILE='diff_tabnew.param.txt.sht', FORM='formatted',
     +      STATUS= 'OLD',IOSTAT=IRR,ERR=60 )
            write(6,*) ' read partons from file: diff_tabnew.param.txt'
         ELSEIF(NG.EQ.16) THEN
            open(30,FILE='h1flat.qg.dat', FORM='unformatted',STATUS=
     +      'OLD',IOSTAT=IRR,ERR=60 )
            write(6,*) ' read partons from file: h1flat.qg.dat'
         ELSEIF(NG.EQ.17) THEN
            open(30,FILE='h1flat.q.dat', FORM='unformatted',STATUS=
     +      'OLD',IOSTAT=IRR,ERR=60 )
            write(6,*) ' read partons from file: h1flat.q.dat'
         ELSE
            write(6,*) ' requested parton distribution not available'
            write(6,*) ' PROGRAM STOP !!!!!!!!!!!!!!!!!!!!!'
            STOP
         ENDIF
         READ(30,'(A)') Line
*         READ(30) XMIN,XMAX,Q2MIN,Q2MAX,NX,NQ

*         write(6,*) 'XMIN,XMAX,Q2MIN,Q2MAX,NX,NQ',
*     +    XMIN,XMAX,Q2MIN,Q2MAX,NX,NQ
         NX=100
cold     NQ=30
ctabnew     NQ=80
         NQ = 50
c xnorm comes from my norm of streng flux =58.4/16/pi=1.169
         xnorm = 1.D0/1.169D0
         DO 20 J=1,NQ
            DO 10 I=0,NX-1
cold               READ(30,*) XX(I),Q2X(J),XGLU,XPART
               READ(30,*) XX(I),Q2X(J),XPART,XGLU,xcharm
c  c.s.f is F_2_charm...
c         write(6,*) ' XX,Q2X,xGlu,xpart',xx(i),q2x(j),xglu,xpart
               do k=-3,3
                  XPD(I,J,K) = xpart/6.D0*xnorm
                  if(k.eq.0) XPD(I,J,K) = xglu*xnorm
               enddo
               XPD(I,J,-4) = xcharm*xnorm/2.D0/9.D0/4.D0
               XPD(I,J,4) = XPD(I,J,-4)
c         write(6,*) ' XX,Q2X,XPD',I,J,
c     &   XX(I),Q2X(I),(XPD(I,J,K),K=-6,6)
   10       CONTINUE
            XX(100)=1.D0
            do k=-6,6
               XPD(I,J,K) = 0.D0
            enddo
   20    Continue

c         write(6,*) (xx(i),i=0,nx)
c         write(6,*) (q2x(i),i=1,nq)
c         do 22 i=0,nx
c         do 22 j=1,nq
c  22     write(6,*) (xpd(i,j,k),k=-6,6)
         FIRST=.FALSE.
         write(6,*) ' laurent pomeron pdf read from file unit 30 '
      ENDIF


      XPRT = DBLE(X)
      Q2T= DBLE(Q2)
      IF(Q2T.LT.Q2X(1)) Q2T=Q2X(1)
      IF(XPRT.LT.XX(0).OR.XPRT.GT.XX(NX) .OR.Q2T.LT.Q2X(1) .OR.Q2T
     +.GT.Q2X(NQ)) THEN
         IF(npr.lt.5) THEN
            WRITE(6,*) 'RES : X or Q2 values outside grid '
            WRITE(6,*) ' X_min ',XX(0),' X_max ',XX(NX), ' actual '
     +      //'X ', XPRT
            WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NQ), ' ' //
     +      'actual Q2 ',Q2T
         ENDIF
         IF(XPRT.LT.XX(0)) XPRT=XX(0)
         IF(XPRT.GT.XX(NX)) XPRT=XX(NX)
         IF(Q2T.LT.Q2X(1)) Q2T = Q2X(1)
         IF(Q2T.GT.Q2X(NQ)) Q2T = Q2X(NQ)
         npr = npr + 1
      ENDIF
      IX = -1
   30 IX = IX + 1
      IF(IX.GT.NX-1) write(6,*) IX,XPRT
      IF(XPRT.GT.XX(IX+1)) GOTO 30
      IQ =  0
   40 IQ = IQ + 1
      IF(IQ.GT.NQ-1) write(6,*) IQ,Q2T
      IF(Q2T.GT.Q2X(IQ+1)) GOTO 40
      DO 50 IP=-6,6
         XPQ(IP)=0.0
         XD = (XPRT - XX(IX))/(XX(IX+1)-XX(IX))
         QD = (Q2T - Q2X(IQ))/(Q2X(IQ+1) - Q2X(IQ))
         X1P=(XPD(IX+1,IQ,IP)-XPD(IX,IQ,IP))*XD +XPD(IX,
     +   IQ,IP)
         X2P=(XPD(IX+1,IQ+1,IP)-XPD(IX,IQ+1,IP))*XD + XPD(
     +   IX,IQ+1,IP)
         XPQ(IP) = SNGL((X2P-X1P)*QD + X1P)
c         write(6,*) x2p,x1p,qd
   50 CONTINUE
c         write(6,*) 'x,q2,xpq',x,q2,(xpq(i),i=-3,3)
      return
   60 write(6,*) ' error in opening file '
      stop
      END
