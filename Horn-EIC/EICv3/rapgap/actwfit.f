*CMZ :  2.07/03 26/03/99  17.10.39  by  Hannes Jung
*CMZ :  2.07/02 24/03/99  16.25.34  by  Hannes Jung
*CMZ :  2.07/01 17/03/99  10.54.12  by  Hannes Jung
*-- Author :
      SUBROUTINE ACTWFIT(BETA_IN,Q2_IN,XPQ,X_POM_IN,T2_IN)

      IMPLICIT Double Precision (A-G,O-Z)
      REAL*4 T2_IN,X_POM_IN,BETA_IN,Q2_IN,XPQ(-6:6)
      Double Precision T2,X_POM,DELTA2,MP,ALPHAP,ALPHA0
      DATA DELTA2/3.24D0/,MP/0.938D0/,ALPHAP/0.25D0/
      LOGICAL FIRST
      INTEGER ICALL
      DATA ICALL/0/
      LOGICAL FLAG
      DATA FLAG/.FALSE./
      integer iret, i
      character*(78) header, filename
      Parameter (Isetmax=15)
      character Flnm(Isetmax)*13
      Data (Flnm(I), I=1,Isetmax)
     + /'ACTW_A2.dat  ','ACTW_B2.dat  '
     + ,'ACTW_C2.dat  ','ACTW_D2.dat  ','ACTW_SG2.dat '
     + ,'ACTW_A2+.dat ','ACTW_B2+.dat '
     + ,'ACTW_C2+.dat ','ACTW_D2+.dat ','ACTW_SG2+.dat'
     + ,'ACTW_A2-.dat ','ACTW_B2-.dat '
     + ,'ACTW_C2-.dat ','ACTW_D2-.dat ','ACTW_SG2-.dat'/
      double precision pardis
*KEEP,RGDIFFR.
      INTEGER NG,NPOM
      DOUBLE PRECISION T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA
      COMMON/DIFFR/T2MAX,XF,ALPHP,RN2,EPSP,QMI,YMI,QMA,YMA,NG,NPOM
      INTEGER IREM
      COMMON/PREMNANT/IREM
*KEND.

      DATA FIRST/.TRUE./
      Icall = icall + 1
      if(icall.eq.1) then
         first=.true.
      else
         first=.false.
      endif

C Inform user of folly
      IF(FIRST) THEN
         WRITE(6,*)'#############################################'
         WRITE(6,*)'#            ACTW fit Selected              #'
         if(ng.eq.-1) then
            write(6,*)'#            fit A                          #'
         elseif(ng.eq.-2) then
            write(6,*)'#            fit B                          #'
         elseif(ng.eq.-3) then
            write(6,*)'#            fit C                          #'
         elseif(ng.eq.-4) then
            write(6,*)'#            fit D                          #'
         elseif(ng.eq.-5) then
            write(6,*)'#            fit SG (singular gluon)        #'
         else
            write(6,*)'#  no valid fit selected (-5 < NG < -1)     #'
            write(6,*)'#            program stopped                #'
            STOP
         endif
         IF(NPOM.EQ.-1) THEN
            write(6,*)'#         DL flux with alpha(0)=1.144       #'
         ELSEIF(NPOM.EQ.-2) THEN
            write(6,*)'#         DL flux with alpha(0)=1.189       #'
         ELSEIF(NPOM.EQ.-3) THEN
            write(6,*)'#         DL flux with alpha(0)=1.085       #'
         ELSE
            write(6,*)'#  no valid fit selected (-3 < NPOM < -1)   #'
            write(6,*)'#            program stopped                #'
            STOP
         endif
         WRITE(6,*)'#############################################'
         IF(NPOM.EQ.-1) THEN
            ALPHA0 = 1.144D0
         ELSEIF(NPOM.EQ.-2) THEN
            ALPHA0 = 1.189D0
         ELSEIF(NPOM.EQ.-3) THEN
            ALPHA0 = 1.085D0
         ELSE
            write(6,*) ' npom = ',npom,' not implemented '
            write(6,*) ' valid values are npom = -1,-2,-3'
            write(6,*) ' program stooped now '
            stop
         endif
         Iset = ABS(NG) + (ABS(NPOM)-1)*5
cwrong      call SetCtq4 (Iset)
         filename=Flnm(Iset)
         Print* ,'Opening ', filename
         call evlrnfe(filename, header, iret)


      ENDIF

C Input quantities are REAL*4 -> Convert to REAL*8
      BETA =DBLE(BETA_IN)
      Q    =DBLE(Sqrt(Q2_IN))
      X_POM=DBLE(X_POM_IN)
      T2   =DBLE(T2_IN)

      PI = 4.D0*DATAN(1.D0)
      F_POM = 9.D0*DELTA2/4.D0/PI**2
      F_DL = (4.D0*MP**2 - 2.8D0*T2)/(4.D0*MP**2 - T2)
      F_DL = F_DL/(1.D0 - T2/0.7D0)**2
      F_DL = F_DL**2
      F_POM = F_POM*F_DL*(X_POM**(1.D0-2.D0*(ALPHA0 + ALPHP*T2)))

      DO I=-6,6
         XPQ(I)=0.
cccc         XPART = Real(Beta*Ctq4Pdf (I, Beta, Q))*F_POM
         XPART = Real(Beta*pardis(I, Beta, Q, iret))*F_POM
         IF(XPART.EQ.0D0) XPART=1.D-10
         XPQ(I) = XPART
cc      write(6,*) ' I ',I,' xpdf ',Beta*pardis(I, Beta, Q, iret),F_POM
      ENDDO

      FIRST=.FALSE.

      RETURN
      END
