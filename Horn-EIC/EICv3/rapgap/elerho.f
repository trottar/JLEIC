*CMZ :  2.08/04 22/12/99  15.39.29  by  Hannes Jung
*CMZ :  2.06/34 22/05/98  17.51.25  by  Hannes Jung
*CMZ :  2.06/02 29/07/97  09.17.51  by  Hannes Jung
*CMZ :  2.03/02 16/08/96  10.33.55  by  Hannes Jung
*CMZ :  2.01/09 05/03/96  09.53.48  by  Hannes Jung
*CMZ :  2.01/04 24/01/96  16.10.48  by  Hannes Jung
*CMZ :  2.01/03 20/01/96  17.12.39  by  Hannes Jung
*CMZ :  2.01/01 09/01/96  14.38.40  by  Hannes Jung
*-- Author :

      SUBROUTINE ELERHO(WT1)
C
C    POMERON-VIRTUAL GAMMA ----> POMERON*+RHO
C    double diffractive ineraction
C
      IMPLICIT REAL*8 (A-H,O-Z)
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

*KEEP,RGLUCO.
      REAL PLEPIN,PPIN
      INTEGER KE,KP,KEB,KPH,KGL,KPA,NFRAG,ILEPTO,IFPS,IHF,IALMKT
      INTEGER INTER,ISEMIH
      INTEGER NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT,NFLAV,NFLQCDC
      COMMON/LUCO  /KE,KP,KEB,KPH,KGL,KPA,NFLAV,NFLQCDC
      COMMON/INPU  /PLEPIN,PPIN,NFRAG,ILEPTO,IFPS,IHF,IALMKT,INTER,
     +              ISEMIH
      COMMON/HARD/ NIA1,NIR1,NIA2,NIR2,NF1,NF2,NFT
      INTEGER IHFLA
      COMMON/HFLAV/ IHFLA
C      SAVE

*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEND.

      DIMENSION PPG(10,5),BETA1(3),BETA2(3),
     +POMSDC(3,5)
      DOUBLE PRECISION P0(4),PG1(4),PG2(4)
      DOUBLE PRECISION SPHI,STHETA,THETA1,PHI1,THETA2,PHI2,
     +POMSPG(3),pp(2,5),
     +T46PG,XM1L,XM1H,Y,STETAS,TETA2,
     +CTETA2,XM1,XM2,XM1S,XM2S,POMNOR
      REAL XMPOMS,XMPIP,XMPIM,XMPI0,POMSDC
      COMMON/DECA2/PG1,PG2
      LOGICAL FIRST
      Double Precision draprn
      DATA FIRST/.TRUE./
      DATA NEV,IEVENT,NEVNT,NSTOP/0,1,10,200/
      DATA NEVE/0/
      DATA XM1L,XM1H/3.000,3.100/
      DATA XMPIP,XMPIM,XMPI0/0.13956755,0.13956755,0.1349734/
      XMPIS=0.0194791D0
      IF(FIRST) THEN
C     write(6,*) ' WE ARE IN ELERHO gamma* pomeron*--> rho pomeron '
         FIRST = .FALSE.
         DO J=1,10
            DO I=1,5
               PPG(J,I) = 0.
            ENDDO
         ENDDO
      ENDIF
      NEVE=0
C...all of the following vectors are or will be given in the MAIN PG cms
C...  P(1, ) is the incoming electron
c...  p(2, ) the incoming proton
c...  p(3, ) the scattered electron
c...  p(nia1, ) the virtual gamma
c...  p(nia2, ) the POMERON
c...  with nia1=4 and nia2=5
c...  p(6, ) is the RHO
C...  P(7, ) is the  EXITED POMERON
c...  p(8, ) is the scattered proton
      DO J=1,5
         P(8,J) = DBLE(P(2,J)-P(NIA2,J))
      ENDDO
      IF(IEVENT.GT.3) GO TO 10
c      PRINT *,' P(1,1)=',P(1,1),' P(1,2)=',P(1,2),' P(1,3)=',P(1,3)
c      PRINT *,' P(2,1)=',P(2,1),' P(2,2)=',P(2,2),' P(2,3)=',P(2,3)
c      PRINT *,' P(3,1)=',P(3,1),' P(3,2)=',P(3,2),' P(3,3)=',P(3,3)
c... the P(NIA2,J) being the pomeron
c      write(6,*)' elerho: nia1,nia2',nia1,nia2,' N=',n
c      PRINT *,' original POMERON end GAMMA in GP cms'
c      PRINT *,' P(5,1)=',P(5,1),' P(5,2)=',P(5,2),' P(5,3)=',P(5,3)
c      PRINT *,' P(4,1)=',P(4,1),' P(4,2)=',P(4,2),' P(4,3)=',P(4,3)
   10 CONTINUE
c
c... so,here the vectors P(NIA1,J) and P(NIA2,J) are given in the PG cms
c...but the pomeron is not  aligned with the z axis.We will now turn it, but
c...before that we conserve their values:
      DO I=1,5
         PP(1,I)=P(NIA1,I)
         PP(2,I)=P(NIA2,I)
      ENDDO
c
c
c notice the combined order in wich ULANGL and LUDBRB must be used !!!!!!
c
      SPHI=DLANGL(P(NIA2,1),P(NIA2,2))
      CALL DUDBRB(NIA1,NIA2,0.D0,-SPHI,0.D0,0.D0,0.D0)
      STHETA=DLANGL(P(NIA2,3),P(NIA2,1))
      CALL DUDBRB(NIA1,NIA2,-STHETA,0.D0,0.D0,0.D0,0.D0)
      IF(IEVENT.LE.0) THEN
         PRINT *,' CM(1)=',CM(1),' CM(2)=',CM(2),' CM(3)=',CM(3)
         PRINT *,' CM(4)=',CM(4)
         PRINT *,' DBCMS(1)=',DBCMS(1),' DBCMS(2)=',DBCMS(2)
         PRINT *,' DBCMS(3)=',DBCMS(3),' DBCMS(4)=',DBCMS(4)
      ENDIF
      IF(IEVENT.GT.0) GO TO 20
      PRINT *,' THE INITIAL ANGLE OF THE POMERON IN PG cms'
      PRINT *,' STHETA=',STHETA,' SPHI=',SPHI
      PRINT *,' APRES ROTATION'
      PRINT *,' P(5,1)=',P(5,1),' P(5,2)=',P(5,2),' P(5,3)=',P(5,3)
      PRINT *,' P(4,1)=',P(4,1),' P(4,2)=',P(4,2),' P(4,3)=',P(4,3)
c      CALL DULIST(1)

   20 CONTINUE
c
C
c...... remember the angle STHETA and SPHI
c*********************************************************************
c
C
C... IN THE POMERON-VIRTUAL GAMMA (MAIN PG CMS) CMS
C... THE VECTORS P(5, ) ,POMERON, AND P(4, ), VIRTUAL GAMMA, ARE
C... DEFINED WITH RESPECT TO AXIS PARALELES TO THE EP CMS SYSTEM BUT TRAVELING
C... WITH A VELOCITY VECTOR GIVEN BY :  -dbcms(1)/dbcms(4),-dbcms(2)/dbcms(4)
C...  -dbcms(3)/dbcms(4) . )do not forget the - sign !!!!!!!!!!)
C..  IN THIS MAIN PG CMS , THE POMERON IS NOW ALIGNED WITH THE Z AXIS.
c
c
C..  THE DOUBLE-DIFFRACTIVE INTERACTION
C...     POMERON-VIRTUAL GAMMA----> POMERONSTAR + RHO
C... FOLLOWED BY THE DECAY TO two pions OF THE POMERON*  AND THE RHO
C...-----------------------------------------------------------------
C
C
C... RELATIVE TO THE POMERON'S DIRECTION OF MOTION (which is actualy  the
c  z axis)
C... THE DIFFERENCIAL CROSS-SECTION sdd(S)/(dT*dXM1S*dXM2S) IS
C... WRITTEN AS = A*1./(XM1S*XM2S)*EXP(B*T)*FDD, if we use non linear mass axis
C...or A*1./(XM1*XM2)*EXP(B*T)FDD, if we use linear mass axis i.e.
C....we write the cross-section as:sdd(S)/(dT*dXM1*dXM2) using dXM1 in place of
C...dXM1S and dXM2 in place of dXM2S
C... WHERE T : IS THE SQUARED MOMENTUM TRANSFER BETWEEN POM AND POM*
c... and the T will be used to obtain the angular distribution of POM* in the
c. . MAIN PG CMS REFERENCE.
C...      XM1S: IS THE SQUARED MASS OF POMERON STAR
C...      XM2S: IS THE SQUARED MASS OF GAMMA STAR
C...      A : CONSTAN
C... FOR OUR PURPOSE , WE CONSIDER FDD ALSO AS CONSTANT . THE MASSES
C... XM1 AND XM2 CAN BE HANDLED AS FOLLOWS : FOR XM1 (THE MASS OF THE
C... EXCITED POMERON) WE WILL CHOSE values starting from 0.420 GeV (3*pion)
c... and following a 1/XM1 pattern .(i.e.we use linear mass axis)
C... FOR XM2 (THE MASS OF RHO ) WE WILL USE THE BREIT-WIGNER . OUR AIM
C... BEING TO GENERATE THE RIGHT T DISTRIBUTION ( FOR GIVEN VALUES OF
C... XM1 AND XM2) , WE CONCENRATE  ONLY ON EXP(BDD*T) and hand xm1 also nearly
C... fixed.





C...choice of the masses-----
      SHAT=(P(4,4)+P(5,4))**2
c      SHAT1=P(4,4)**2-(P(4,1)**2+P(4,2)**2+P(4,3)**2)+
c     &P(5,4)**2-(P(5,1)**2+P(5,2)**2+P(5,3)**2)+2.*P(4,4)*P(5,4)-
c     &2.*(P(4,1)*P(5,1)+P(4,2)*P(5,2)+P(4,3)*P(5,3))
c      IF(SHAT.NE.SHAT1) PRINT *,' SHAT=',SHAT,' SHAT1=',SHAT1
      ETOT=SQRT(SHAT)
C     PRINT *,' SHAT=',SHAT
c
      CALL HFILL(606,SNGL(ETOT),0.,1.)
C
C... for the rho
C
   30 X=draprn()
      NEVE=NEVE+1
      IF(NEVE.GE.NSTOP) WT1=0.D0
c      IF(NEVE.GE.NSTOP) PRINT *,' BBBBBBB'
      IF(NEVE.GE.NSTOP) CALL HFILL(932,1,0.,1.)
      IF(NEVE.GE.NSTOP) GO TO 70
      XM2=0.770+X*0.075D0
      Y=draprn()
      IF(Y.LT.0.5) XM2=0.770-X*0.075D0
      XM2S=XM2*XM2
      IF(XM2.LT.0.620.OR.XM2.GT.0.920) GO TO 30
C     PMAX=SQRT((SHAT-XM1MIN**2-XM2**2)**2/4.-(XM1MIN*XM2)**2)/ETOT
C     PRINT *,' XM2=',XM2

      CALL HFILL(605,SNGL(XM2),0.,1.)
C
c...for the pomeron*
c
c
c....if XM1 is not fixed
c 112  Y=draprn()
c      XM1=(1./Y)*0.420
c      IF(XM1.GT.3.) GO TO 112
c      CALL HFILL(963,SNGL(XM1),0.,1.)


c...if XM1 is fixed
c
   40 Y=draprn()
      IF(Y.EQ.0.) GO TO 40
      NEVE=NEVE+1
      IF(NEVE.GE.NSTOP) WT1=0.D0
c       IF(NEVE.GE.NSTOP) PRINT *,' AYAYA'
      IF(NEVE.GE.NSTOP) CALL HFILL(932,1,0.,1.)
      IF(NEVE.GE.NSTOP) GO TO 70
      XM1=XM1L+Y*(XM1H-XM1L)
      IF(XM1.GE.ETOT) GO TO 40
C     PRINT *,' CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
C     PRINT *,' CCCCCCCCCCCCCCCC XM1=',XM1,' CCCCCCCCC'
C     PRINT *,' CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      XM1S=XM1*XM1
c
c
      XTOT=XM1+XM2
      IF(XTOT.GT.ETOT) GO TO 40
      PHAT=(SHAT+XM1S-XM2S)**2/(4.D0*SHAT)-XM1S
      if(PHAT.LT.0.) GO TO 30
      PHAT=SQRT(PHAT)
C...generate the interaction T (impulsion transfer) distribution.
C
C... the choice of the excited masses imposes
C... limits on the inpusion transfer, T ,
C...
C
C     PRINT *,' SHAT=',SHAT
      XMU1=-PPG(5,5)**2/SHAT
      XMU2=-PPG(6,5)**2/SHAT
      XMU3=XM1S/SHAT
      XMU4=XM2S/SHAT
C     PRINT *,' XMU1=',XMU1,' XMU2=',XMU2
C     PRINT *,' XMU3=',XMU3,' XMU4=',XMU4
      C1=1.D0-(XMU1+XMU2+XMU3+XMU4)+(XMU1-XMU2)*(XMU3-XMU4)
      C111=(1.D0-XMU1-XMU2)**2-4.*XMU1*XMU2
      C222=(1.D0-XMU3-XMU4)**2-4.*XMU3*XMU4
      IF(C111.LT.0.OR.C222.LT.0.) PRINT 10000
10000 FORMAT(' PROBLEMS WITH C1 OR C2')
      C2=SQRT(C111)*SQRT(C222)
      C3=(XMU3-XMU1)*(XMU4-XMU2)+(XMU1+XMU4-XMU2-XMU3)*(XMU1*XMU4-
     +XMU2*XMU3)
C     PRINT *,' C1=',C1,' C2=',C2
C     PRINT *,' C3=',C3
      TMIN=-SHAT*(C1+C2)/2.D0
      TMAX=-SHAT*(C1-C2)/2.D0
C     PRINT *,' TMIN=',TMIN,' TMAX=',TMAX
C
C... we must have TMIN<T<TMAX
C

      BDD=0.5*LOG(2.72D0**4+4.D0*SHAT/(XM1S*XM2S))
      YY=BDD*TMIN
C     PRINT *,' YY=',YY
      IF(YY.GT.-78.) THEN
         YMIN=EXP(BDD*TMIN)
         YMAX=EXP(BDD*TMAX)
C     PRINT *,' YMIN=',YMIN,' YMAX=',YMAX
         Y=YMIN+(YMAX-YMIN)*draprn()
C     PRINT *,' Y=',Y,' BDD=',BDD
         T=LOG(Y)/BDD
         GO TO 50
      ENDIF
      YMIN=0.D0
      YMAX=EXP(BDD*TMAX)
C     PRINT *,' YMIN=',YMIN,' YMAX=',YMAX
      Y=YMIN+(YMAX-YMIN)*draprn()
C     PRINT *,' Y=',Y,' BDD=',BDD
      T=LOG(Y)/BDD
   50 CONTINUE
c...****************************************************
c.. energy of the POMERON* and of the RHO
C...****************************************************
      P(7,4)=SQRT(PHAT**2+XM1S)
      P(6,4)=SQRT(PHAT**2+XM2S)
      EPRIE1=P(7,4)
      EPRIE2=P(6,4)
C*******************************************************
C
C...
C... using the the generated T value , get the cosinus of the the angle theta
C ..  of the pomeron* (in the MAIN PG cms
c...
      POMNOR=SQRT(P(5,1)**2+P(5,2)**2+P(5,3)**2)
      COSTH=(P(5,4)*P(7,4)+(T+P(5,5)**2-XM1S)/2.D0)/
     +(POMNOR*PHAT)
c      IF(COSTH.LT.0.) THEN
c      PRINT *,' COSTH NEGATIVE'
c      PRINT *,' XM1=',XM1,' P(5,5)=',P(5,5),' T=',T,' PHAT=',PHAT
c      ENDIF
      IF(COSTH.GE.1.) GO TO 30
      CALL HFILL(609,SNGL(EPRIE1),0.,1.)
      CALL HFILL(610,SNGL(EPRIE2),0.,1.)
      CALL HFILL(630,SNGL(PHAT),0.,1.)
      CALL HFILL(607,SNGL(XM1),0.,1.)
      CALL HFILL(608,SNGL(XTOT),0.,1.)
      CALL HFILL(631,SNGL(T),0.,1.)

C     PRINT *,' T=',T,' EPRIE1=',EPRIE1,' EPRIE2=',EPRIE2
c      PRINT *,' POMNOR=',POMNOR,' PHAT=',PHAT
c      PRINT *,' P(5,4)=',P(5,4),' P(5,5)=',P(5,5),
c    &' P(7,4)=',P(7,4)
c      PRINT *,' T=',T,' XM1S=',XM1S
      XTEST=P(5,5)
      CALL HFILL(800,SNGL(COSTH),POMNOR,1.)
      CALL HFILL(801,SNGL(COSTH),SNGL(PHAT),1.)
      CALL HFILL(802,SNGL(COSTH),SNGL(XTEST),1.)
      CALL HFILL(803,SNGL(COSTH),SNGL(T),1.)
      CALL HFILL(680,SNGL(ETOT),0.,1.)
C     PRINT *,' T=',T,' COSTH=',COSTH
      CALL HFILL(632,SNGL(COSTH),0.,1.)
C
C
C
C... excited POMERON's angles in the MAIN PG cms
***********************************************************************
C************** remember these angles !!!!!!!!!!!!!!!!!!!!!!!!!!!
      THETA1=ACOS(COSTH)
      PHI1=draprn()*6.28D0
      IF(IEVENT.LE.0)PRINT *,'POM* ANGLES FROM GENERATOR IN MAIN PG cms'
      IF(IEVENT.LE.0)PRINT *,' THETA1=',THETA1,' PHI1=',PHI1
c...
      CALL HFILL(633,SNGL(THETA1),0.,1.)
c... the POM*  four momentum in the reference MAIN PG cms
c..
c..
      P(7,1)=PHAT*SIN(THETA1)*COS(PHI1)
      P(7,2)=PHAT*SIN(THETA1)*SIN(PHI1)
      P(7,3)=PHAT*COS(THETA1)
C...FOR LATTER USE
      POMSPG(1)=P(7,1)
      POMSPG(2)=P(7,2)
      POMSPG(3)=P(7,3)
C
c...  p(7,4) is already defined
      P(7,5)=SQRT(P(7,4)**2-(P(7,1)**2+P(7,2)**2+P(7,3)**2))
c...     the RHO's  four momentum in MAIN PM cms
      P(6,1)=-P(7,1)
      P(6,2)=-P(7,2)
      P(6,3)=-P(7,3)
C      P(6,4) already defined
      P(6,5)=SQRT(P(6,4)**2-(P(6,1)**2+P(6,2)**2+ P(6,3)**2))
c..    ......the velocity vector of the POM* in the MAIN PG cms is now
      BETA1(1)=P(7,1)/P(7,4)
      BETA1(2)=P(7,2)/P(7,4)
      BETA1(3)=P(7,3)/P(7,4)
C
C... the velocity vector of the RHO in the MAIN PG cms
C
      BETA2(1)=P(6,1)/P(6,4)
      BETA2(2)=P(6,2)/P(6,4)
      BETA2(3)=P(6,3)/P(6,4)
c...the gamma*-rho momentum transfer,in quantities expressed in PG cms;
      T46PG=P(4,4)**2-(P(4,1)**2+P(4,2)**2+P(4,3)**2)+XM2S-2.D0*P(4,4)*
     +P(6,4)+2.D0*(P(4,1)*P(6,1)+P(4,2)*P(6,2)+P(4,3)*P(6,3))
      CALL HFILL(927,SNGL(T46PG),0.,1.)
      IF(IEVENT.LE.0) THEN
         PRINT *,' the velocity vector of pomeron* in the MAIN PG cms'
         PRINT *,' BETA1(1)=',BETA1(1),BETA1(2),BETA1(3)
         PRINT *,' the velocity vector of rho  in the MAIN PG cms'
         PRINT *,' BETA2(1)=',BETA2(1),BETA2(2),BETA2(3)
      ENDIF

c...in the following NF1 corresponds to the RHO and NF2 to POM*
C....in our case NF1=6 end NF2=7
C
      IF(IEVENT.LE.0) PRINT *,' NF1=',NF1,' NF2=',NF2
c
c
c.... DECAY OF THE POM* AND RHO
C....******************************
c    in the POM* cms
C... we rotate  the pom* flighing path  to the z axis
c...
c
      IF(IEVENT.LE.0)PRINT * ,' original POM* end RHO in PG cms'
      IF(IEVENT.LE.0)PRINT *,' p(7,1)=',p(7,1),' p(7,2)=',p(7,2),
     +' p(7,3)=',p(7,3)
      IF(IEVENT.LE.0)PRINT *,' p(6,1)=',p(6,1),' p(6,2)=',p(6,2),
     +' p(6,3)=',p(6,3)
      PHI1=DLANGL(P(7,1),P(7,2))
      CALL DUDBRB(7,7,0.D0,-PHI1,0.D0,0.D0,0.D0)
      THETA1=DLANGL(P(7,3),P(7,1))
      CALL DUDBRB(7,7,-THETA1,0.D0,0.D0,0.D0,0.D0)
      IF(IEVENT.LE.0) THEN
         PRINT *,' THETA1 AND PHI1 FROM ULANGL'
         PRINT *,' THETA1=',THETA1,' PHI1=',PHI1
      ENDIF
      IF(IEVENT.LE.0)PRINT * ,' POM* and RHO IMP AFTER ROTATION'
      IF(IEVENT.LE.0)PRINT *,' P(7,1)=',P(7,1),' P(7,2)=',P(7,2),
     +' P(7,3)=',P(7,3)
c...remember that for the reverse rotation we must first rotate THETA1
C...followed by the rotation PHI1
C
C..
C.. we do the same for rho (in the rho cms
c ..
c
      PHI2=DLANGL(P(6,1),P(6,2))
      CALL DUDBRB(6,6,0.D0,-PHI2,0.D0,0.D0,0.D0)
      THETA2=DLANGL(P(6,3),P(6,1))
      CALL DUDBRB(6,6,-THETA2,0.D0,0.D0,0.D0,0.D0)
      IF(IEVENT.LE.0)PRINT *,' P(6,1)=',P(6,1),' P(6,2)=',P(6,2),
     +' P(6,3)=',P(6,3)
      IF(IEVENT.LE.0)PRINT *,' ANGLES OF RHO'
      IF(IEVENT.LE.0)PRINT *,' THETA2=',THETA2,' PHI2=',PHI2
c
c
C
C...
C
      K(NF1,2) = 113
      K(NF1,1) = 11
      K(NF2,2) = 100
      K(NF2,1) = 11
c      IF(IEVENT.eq.9)  then
c      PRINT *,' BEFOR THE JUMP'
c      stop 9
c      endif
C
C... the RHO being polarised
C...  and therefor its decay angular distribution
C... in RHO giving PI+PI- is taken to be proportional to SIN(TETA)**2
C... where the reference axis is given by the rho direction of motion.
C... in this case the z axis
c   ******************************************************************
c ... For the POM* we supose no such polarisation.
c....
c... *****************************************************************
c
c..... decay of the POM* to three pions (pi+,pi-,p0)
c
      IF(IEVENT.LE.0) PRINT *,' BEGINS DECAY GDECA3 IN THREE PIONS'
      IF(IEVENT.LE.0) PRINT *,' P(7,5)=',P(7,5)
      XMPOMS=SNGL(P(7,5))
c hju line commented because XM1MIN was not defined
c      IF(XMPOMS.LT.XM1MIN) PRINT *,' PROBLEME WITH POM* MASS'


c
c
      CALL GDECA3(XMPOMS,XMPIP,XMPIM,XMPI0,POMSDC)
      IF(IEVENT.LE.0) PRINT *,' HAVE PASSED GDECA3'
      IF(IEVENT.LE.0) PRINT *,' POMSDC(1,1)=',POMSDC(1,1),
     +' POMSDC(1,2)=',POMSDC(1,2),' POMSDC(1,3)=',POMSDC(1,3)
C
C... the pion plus in the pom* cms
C
      NPI = N + 1
      K(NPI,1)=1
      K(NPI,2)=211
      K(NPI,3)=7
      DO I=1,5
         P(NPI,I)=DBLE(POMSDC(1,I))
      ENDDO
      IF(IEVENT.LE.0) THEN
         PRINT *,' the pion plus NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3),' P4='
     +   , P(NPI,4)
      ENDIF
c
c.... the pi- in the POM* cms
c
      NPI = NPI + 1
      K(NPI,1)=1
      K(NPI,2)=-211
      K(NPI,3)=7
      DO I=1,5
         P(NPI,I)=DBLE(POMSDC(2,I))
      ENDDO
      IF(IEVENT.LE.0) THEN
         PRINT *,' the pion minus NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
c
c....the pi0 in the POM* cms
c
      NPI = NPI + 1
      K(NPI,1)=11
      K(NPI,2)=111
      K(NPI,3)=7
      NPI0 = NPI
      DO I=1,5
         P(NPI,I)=DBLE(POMSDC(3,I))
      ENDDO
      IF(IEVENT.LE.0) THEN
         PRINT *,' the pi zero NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF


c... and after that,rotate them back with the angles of the POM*
c... in the MAIN PG cms
      CALL DUDBRB(NPI-2,NPI,THETA1,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(NPI-2,NPI,0.D0,PHI1,0.D0,0.D0,0.D0)
c
c... rotate back also the POM* vector and check that we get the old value
c
      CALL DUDBRB(7,7,THETA1,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(7,7,0.D0,PHI1,0.D0,0.D0,0.D0)

C
C...
      IF(IEVENT.LE.0) THEN
         PRINT *,'  AFTER BACKWARD ROTATION WITH THE POM* ANGLES '
         PRINT *,' FOR CHECK : THE POMERON'
         PRINT *,' P(7,1)=',P(7,1),' P(7,2)=',P(7,2),' P(7,3)=',P(7,3)
         PRINT *,' THE PI PLUS NPI=',NPI-2
         PRINT *,' P1=',P(NPI-2,1),' P2=',P(NPI-2,2),' P3=',P(NPI-2,3)
         PRINT *,' THE PI MINUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI ZERO NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
c
c ...boost now the pions only back to the MAIN PG cms
C...
      CALL DUDBRB(NPI-2,NPI,0.D0,0.D0,BETA1(1),BETA1(2),BETA1(3))
      IF(IEVENT.LE.0) THEN
         PRINT *,' the pions , once back in the PG cms'
         PRINT *,' THE PI PLUS NPI=',NPI-2
         PRINT *,' P1=',P(NPI-2,1),' P2=',P(NPI-2,2),' P3=',P(NPI-2,3)
         PRINT *,' THE PI MINUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI ZERO NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
C..
C.. end rotate them back with the angles of the POM in the PG cms
      CALL DUDBRB(NPI-2,NPI,STHETA,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(NPI-2,NPI,0.D0,SPHI,0.D0,0.D0,0.D0)
      CALL DUDBRB(7,7,STHETA,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(7,7,0.D0,SPHI,0.D0,0.D0,0.D0)
C...FOR LATTER USE
      POMSPG(1)=P(7,1)
      POMSPG(2)=P(7,2)
      POMSPG(3)=P(7,3)
      IF(IEVENT.LE.0) THEN
         PRINT *,'  once rotated back with the angles of the POM'
         PRINT *,' THE POMERON*'
         PRINT *,' P(7,1)=',P(7,1),' P(7,2)=',P(7,2),' P(7,3)=',P(7,3)
         PRINT *,' THE PI PLUS NPI=',NPI-2
         PRINT *,' P1=',P(NPI-2,1),' P2=',P(NPI-2,2),' P3=',P(NPI-2,3)
         PRINT *,' THE PI MINUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI ZERO NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF


C...
c..


c..
C
C...decay of the pizero to gamma-gamma in the PG cms
c
      DO I=1,4
         P0(I)=DBLE(P(NPI,I))
      ENDDO
      CALL P0TOGG(1,P0)
C
C...the gamma1 in the PG cms
C
      NPI=NPI+1
      K(NPI,1)=1
      K(NPI,2)=22
      K(NPI,3)=NPI0
      DO I=1,4
         P(NPI,I)=DBLE(PG1(I))
      ENDDO
      P(NPI,5)=0.D0
      IF(IEVENT.LE.0) THEN
         PRINT *,' NPI=',NPI
         PRINT *,' momentum of the photon1 in the PG cms'
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
C
C...the gamma2 in the PG cms
c
      NPI=NPI+1
      K(NPI,1)=1
      K(NPI,2)=22
      K(NPI,3)=NPI0
      DO I=1,4
         P(NPI,I)=DBLE(PG2(I))
      ENDDO
      P(NPI,5)=0.D0
      IF(IEVENT.LE.0) THEN
         PRINT *,' NPI=',NPI
         PRINT *,' momentum of the photon2 in the PG cms'
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF

C..
C... decay of the RHO
C
   60 Y=draprn()
      STETAS=1.D0-Y**2
      TETA2=ASIN(SQRT(STETAS))
      Y=draprn()
      IF(Y.GT.STETAS) GO TO 60
      Y=draprn()
      IF(Y.LE.0.5) TETA2=3.14D0-TETA2
      CTETA2=COS(TETA2)
      CALL HFILL(613,SNGL(CTETA2),0.,1.)
      FI2=draprn()*6.28D0
      PHELP=SQRT(XM2S/4.D0-XMPIS)
C
C... the pion plus in the RHO CMS
C
      NPI=NPI+1
      K(NPI,1)=1
      K(NPI,2)=211
      K(NPI,3)=6
      P(NPI,1)=PHELP*SIN(TETA2)*COS(FI2)
      P(NPI,2)=PHELP*SIN(TETA2)*SIN(FI2)
      P(NPI,3)=PHELP*COS(TETA2)
      P(NPI,4)=XM2/2.D0
      P(NPI,5)=0.13956755D0
C
C
C... the pion minus in the rho cms
C
      NPI=NPI+1
      K(NPI,1)=1
      K(NPI,2)=-211
      K(NPI,3)=6
      P(NPI,1)=-P(NPI-1,1)
      P(NPI,2)=-P(NPI-1,2)
      P(NPI,3)=-P(NPI-1,3)
      P(NPI,4)=XM2/2.D0
      P(NPI,5)=0.13956755D0
      IF(IEVENT.LE.0) THEN
         PRINT *,' the pions in the RHO cms'
         PRINT *,' THE PI PLUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI MINUS NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF

c
C
c
      SBETA2=BETA2(1)**2+BETA2(2)**2+BETA2(3)**2
      IF(SBETA2.GE.1.) PRINT *,' PROBLEME BETA2'
      IF(SBETA2.GE.1.) RETURN
C.. rotate them back with the angles of RHO in the MAIN PG cms :
      CALL DUDBRB(NPI-1,NPI,THETA2,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(NPI-1,NPI,0.D0,PHI2,0.D0,0.D0,0.D0)
      CALL DUDBRB(6,6,THETA2,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(6,6,0.D0,PHI2,0.D0,0.D0,0.D0)

      IF(IEVENT.LE.0) THEN
         PRINT *,' the same, rotated back with the RHO PG cms angles'
         PRINT *,' THE RHO'
         PRINT *,' P(6,1)=',P(6,1),' P(6,2)=',P(6,2),' P(6,3)=',P(6,3)
         PRINT *,' THE PI PLUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI MINUS NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF

c
C...boost them back,the pions only, to the PG cms
      CALL DUDBRB(NPI-1,NPI,0.D0,0.D0,BETA2(1),BETA2(2),BETA2(3))
      IF(IEVENT.LE.0) THEN
         PRINT *,' the same, boosted back to the PG cms'
         PRINT *,' THE PI PLUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI MINUS NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
C..
C.. end rotate them back with the angles SPHI,  STHETA, of the POM
c.. in the PG
      CALL DUDBRB(NPI-1,NPI,STHETA,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(NPI-1,NPI,0.D0,SPHI,0.D0,0.D0,0.D0)
      CALL DUDBRB(6,6,STHETA,0.D0,0.D0,0.D0,0.D0)
      CALL DUDBRB(6,6,0.D0,SPHI,0.D0,0.D0,0.D0)
      IF(IEVENT.LE.0) THEN
         PRINT *,' the same, once rotated back with STHETA and SPHI'
         PRINT *,' THE RHO'
         PRINT *,' P(6,1)=',P(6,1),' P(6,2)=',P(6,2),' P(6,3)=',P(6,3)
         PRINT *,' FOR COMPARISON THE POMERON* IN THE SAME PG cms'
         PRINT *,' POMSPG(1)=',POMSPG(1),' POMSPG(2)=',POMSPG(2),
     +   ' POMSPG(3)=',POMSPG(3)
         PRINT *,' THE PI PLUS NPI=',NPI-1
         PRINT *,' P1=',P(NPI-1,1),' P2=',P(NPI-1,2),' P3=',P(NPI-1,3)
         PRINT *,' THE PI MINUS NPI=',NPI
         PRINT *,' P1=',P(NPI,1),' P2=',P(NPI,2),' P3=',P(NPI,3)
      ENDIF
C
C...end finaly set back the original values for P(NIA1, ) and P(NIA2, )
      DO I=1,5
         P(NIA1,I)=PP(1,I)
         P(NIA2,I)=PP(2,I)
      ENDDO
C
C
c..
      N=NPI
c      CALL DULIST(1)
c      PRINT *,' NEVE=',NEVE
      IEVENT=IEVENT+1
c      PRINT *,' IEVENT=',IEVENT
      WT1 = 1.d0
      RETURN
   70 CONTINUE
      WT1=0.D0
c      PRINT *,' FUNY'
      RETURN
      END
