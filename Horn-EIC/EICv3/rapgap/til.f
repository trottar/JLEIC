      DOUBLE PRECISION FUNCTION TIL(I,L)
      implicit none
      DOUBLE PRECISION PK1,PK2,PK3,PL
      COMMON /QQG/ PK1(4),PK2(4),PK3(4),PL(4)
      DOUBLE PRECISION LK1K2,K1K2,LK1,K1,LK2
      INTEGER I,L,NM
      DOUBLE PRECISION TILT,CF,T1,T2,X2LK1K2,X2K1K2,X2LK1,X2K1,X2LK2
      DOUBLE PRECISION XK22
      DOUBLE PRECISION DK
      EXTERNAL DK
      TILT = 0.D0
      DO NM = 1,2
         CF = 1.D0
         IF(NM.EQ.2) CF = -1.D0
         IF(I.EQ.0) THEN
            LK1K2 = 1.D0
            K1K2 = 1.D0
            LK1 = -1.D0
            K1 = 1.D0
         ELSE
            LK1K2 = CF*PL(I) + PK1(I) + PK2(I)
            K1K2 = PK1(I) + PK2(I)
            LK1 = CF*PL(I) - PK1(I)
            K1 = PK1(I)
         ENDIF
         X2LK1K2 = (CF*PL(1) + PK1(1) + PK2(1))**2+ (CF*PL(2) + PK1(2)
     +   + PK2(2))**2
         X2K1K2 = (PK1(1) + PK2(1))**2 + (PK1(2) + PK2(2))**2
         X2LK1 = (CF*PL(1) - PK1(1))**2 + (CF*PL(2) - PK1(2))**2
         X2K1 = PK1(1)**2 + PK1(2)**2
         T1 = LK1K2/DK(X2LK1K2) + K1K2/DK(X2K1K2) + LK1/DK(X2LK1)
     +   - K1/DK(X2K1)
         LK2 = CF*PL(L) + PK2(L)
         X2LK2 = (CF*PL(1) + PK2(1))**2 + (CF*PL(2) + PK2(2))**2
         XK22 = PK2(1)**2 + PK2(2)**2
         T2 = LK2/X2LK2 - PK2(L)/XK22

         TILT = TILT + T1*T2
      ENDDO
      IF(TILT.LT.0) THEN
c      write(6,*) ' TIL I,L ',TILT,I,L
      ENDIF
      TIL = TILT
      RETURN
      END
