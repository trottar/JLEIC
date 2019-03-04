*CMZ :  2.07/03 24/05/99  17.03.31  by  Hannes Jung
*-- Author :    Hannes Jung   26/12/97
      SUBROUTINE pdisdc
	Implicit None
c perform decay of p diss system
c special treatment to avoid delta prod.
      INTEGER N,K
      REAL P,V
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
*KEEP,RGLUDAT1.
      REAL PARU,PARJ
      INTEGER MSTU,MSTJ
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE

*KEND.
      real ym_str,ym_had
      common/pdisc/ ym_str,ym_had
	integer npdc1,npdc2
      common/pdcyn/npdc1,npdc2
	real P18old,ETT,epzall,xmr2,ym_test
	Integer I,J,IPR,IPT1,IPT2,IPT3,IPAR1,IPAR2,IPAR3,IREP
      if(npdc1.ne.0) then
         K(npdc1,1) = 2
         K(npdc2,1) = 1
      endif
      P18old = PARJ(18)
      PARJ(18)=0.0000001
cccc      call lulist(1)
      call luexec
      PARJ(18) = P18old


c look for p dissociation
      IPR  = -999
      IPT1 = 1
      IPT2 = 1
      IPT3 = 1
      IPAR1 = -999
      IPAR2 = -999
      IPAR3 = -999
      ETT = 0.
      IREP = 0
      ym_str=0.
      ym_had=0.
c      write(6,*) ' checking p - diss '

      DO 10  J=1,5
   10 P(N+1,J) = 0.
      epzall = 0.
      Do 30  I=6,N
         IF(K(I,2).EQ.2210.AND.K(I,3).EQ.2) IPR = I
c       IF(K(I,3).EQ.IPR) THEN
c          DO 490 J=1,5
c 490         P(N+2,J) = P(N+2,J) + P(I,J)
c          write(6,*) ' IPR found',I,IPR
c          ENDIF
         IF(K(I,3).EQ.IPR.AND.IPT1.EQ.1) THEN
            IPT1 = IPT1+1
            IPAR1 = I
c          write(6,*) ' IPT1 found',IPAR1
            IF(K(I,2).EQ.2212) then
               ym_had=P(I,5)
               goto 40
            endif
         ENDIF
cc         IF(K(I,3).EQ.IPAR1.AND.IPT2.EQ.1) THEN
cc            IPT2 = IPT2 + 1
cc            IPAR2 = I
ccc          write(6,*) ' IPT2 found',IPAR2
cc         ENDIF
cccccccc         IF(K(I,3).EQ.IPAR2.AND.IPT3.EQ.1) THEN
         IF(K(I,3).EQ.IPAR1.AND.IPT2.EQ.1) THEN
c here we have the string from p diss
cc            IPT3 = IPT3 + 1
cc            IPAR3 = I
            IPT2 = IPT2 + 1
            IPAR2 = I
cc          write(6,*) ' string IPT3 found',I
            ym_str = real(P(I,5))
         ENDIF
ccc         IF(K(I,3).EQ.IPAR3) THEN
         IF(K(I,3).EQ.IPAR2) THEN
c calculate mass of string from decay products
cc         write(6,*) 'decay particles ',i
            DO 20  J=1,5
   20       P(N+1,J) = P(I,J) + P(N+1,J)
         ENDIF
   30 CONTINUE
      XMR2 = P(N+1,4)**2 - P(N+1,1)**2 - P(N+1,2)**2 - P(N+1,3)**2
      if(xmr2.gt.0.0) then
         ym_had = sqrt(real(xmr2))
      else
         ym_had = -99999.
      endif
   40 continue
ccc       write(6,*) ' my(hadr),my(str) ',ym_had,ym_str
ccc       call lulist(1)
cc      if(ym_had.gt.20.and.ym_had.le.40) then
cc      write(6,*) ' my(hadr),my(str) ',ym_had,ym_str
c       write(6,*) ' P(n+1) ',(P(n+1,i),i=1,4)
c       call lulist(1)
cc       endif
      IF(MSTU(24).NE.0) THEN
         if(npdc1.ne.0) then
            DO 50  J=1,5
   50       P(N+2,J) = P(npdc1,J) + P(npdc2,J)
            ym_test = P(N+2,4)**2-P(N+2,3)**2-P(N+2,2)**2-P(N+2,1)**2
            ym_test = sqrt(ym_test)
         endif
         WRITE(6,*) 'MSTU(24)= ',MSTU(24)
         write(6,*) ' my(hadr),my(str) ',ym_had,ym_str
         write(6,*) ' ym_test ',ym_test
         write(6,*) ' P(n+2) ',(P(n+2,i),i=1,4)
cc         CALL LULIST(1)
      ENDIF
      return
      end
