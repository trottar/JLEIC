*CMZ :  2.08/04 21/12/99  11.08.29  by  Hannes Jung
*CMZ :  4.41/00 13/12/94  10.51.18  by  Hannes Jung
*-- Author :
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DJGVAR(ICHNN,X,Y,Q2)
C---
C   TRANSFER VARIABLES FROM HERACLES TO LEPTO OR OTHER ROUTINE
C   FOR FRAGMENTATION
C---
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      PARAMETER (NMXHEP=2000)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     &                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
      COMMON /HSRESC/ SSCH,Q2SCH,W2SCH,XSCH,YSCH
     &               ,SLIX,SLIY,SLIZ,SLIE,SLIM
     &               ,SLFX,SLFY,SLFZ,SLFE,SLFM

C---DECLARATIONS FOR LEPTO
C---->
C   TRANSER OF KINEMATICAL VARIABLES
C---->
      B=X
      RETURN
      END
