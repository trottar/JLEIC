      SUBROUTINE H1RNV(RVEC,LENV)
************************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN (RCARRY)                                        *
*#                                                                     *
*# ENTRIES ARE:                                                        *
*#     SUBROUTINE  H1RNV(RVEC,LENV) VECTOR OF RANDOM NUMBERS           *
*#     SUBROUTINE  H1RNIN(INSEED)   INITIALISE WITH SEED               *
*#     SUBROUTINE  H1RNIV(IVEC)     INITIALISE/RESTART WITH SEED ARRAY *
*#     SUBROUTINE  H1RNSV(IVEC)     SAVE SEED ARRAY IVEC(25)           *
*#                                                                     *
*# IN ADDITION THERE EXISTS:                                           *
*#     FUNCTION    H1RN()          SINGLE RANDOM NUMBER                *
*#                                                                     *
************************************************************************
      LOGICAL NOTYET
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24),ISEEDS(24),IVEC(25)
      PARAMETER (TWOP12=4096.)
      PARAMETER (ITWO24=2**24,ICONS=2147483563)
      DATA NOTYET /.TRUE./
      DATA I24,J24,CARRY/24,10,0./
      DATA TWOM24/1/
*
*      write(6,*) 'in h1rnv '
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = 314159265
         WRITE(6,'(A,I12)') ' H1RN default initialization: ',JSEED
         TWOM24 = 1.
         DO 10 I= 1, 24
            TWOM24 = TWOM24 * 0.5
            K = JSEED/53668
            JSEED = 40014*(JSEED-K*53668) -K*12211
            IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
            ISEEDS(I) = MOD(JSEED,ITWO24)
   10    CONTINUE
         DO 20 I= 1,24
            SEEDS(I) = REAL(ISEEDS(I))*TWOM24
   20    CONTINUE
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      ENDIF

C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 30  IV= 1, LENV
         UNI = SEEDS(I24) - SEEDS(J24) - CARRY
         IF (UNI .LT. 0.) THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = I24 - 1
         IF (I24 .EQ. 0) I24 = 24
         J24 = J24 - 1
         IF (J24 .EQ. 0) J24 = 24
*       code to avoid exact zeroes
         IF(UNI.EQ.0.)THEN
            UNI=SEEDS(I24)*TWOM24
            IF(UNI.EQ.0.)UNI=2.**(-48)
         ENDIF
         RVEC(IV) = UNI
   30 CONTINUE
      RETURN

************************************************************************

C                    Entry to initialize from one integer
      ENTRY H1RNIN(INSEED)
      NOTYET=.FALSE.
      JSEED = INSEED
      WRITE(6,'(A,I12)') ' H1RN initialized from seed: ',JSEED
      TWOM24 = 1.
      DO 40  I= 1, 24
         TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0) JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   40 CONTINUE
      DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
   50 CONTINUE
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .LT. SEEDS(14)) CARRY = TWOM24
      RETURN

************************************************************************

C           Entry to input and float integer seeds from previous run

      ENTRY H1RNIV(IVEC)
      NOTYET=.FALSE.
      TWOM24 = 1.
      DO 60  I= 1, 24
   60 TWOM24 = TWOM24 * 0.5
      WRITE(6,'(A)') ' Full initialization of H1RN with 25 integers:'
      DO 70  I= 1, 24
         SEEDS(I) = REAL(IVEC(I))*TWOM24
   70 CONTINUE
      CARRY = REAL(MOD(IVEC(25),10))*TWOM24
      ISD = IVEC(25)/10
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = ISD
      RETURN

************************************************************************
C                    Entry to output seeds as integers

      ENTRY H1RNSV(IVEC)
      DO 80  I= 1, 24
         IVEC(I) = INT(SEEDS(I)*TWOP12*TWOP12)
   80 CONTINUE
      ICARRY = 0
      IF (CARRY .GT. 0.)  ICARRY = 1
      IVEC(25) = 1000*J24 + 10*I24 + ICARRY
      RETURN
      END
