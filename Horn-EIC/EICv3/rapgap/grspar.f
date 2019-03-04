*CMZ :  2.05/00 03/03/97  17.00.42  by  Hannes Jung
*-- Author :
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*           G R S - LO - VIRTUAL PHOTON PARAMETRIZATIONS          *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE                  *
*                M. GLUECK, E.REYA, M. STRATMANN :                *
*                    PHYS. REV. D51 (1995) 3220                   *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE EVOLVED PARTONS FOR    *
*        Q**2 / GEV**2  BETWEEN   0.6   AND  5.E4                 *
*                       AND (!)  Q**2 > 5 P**2                    *
*        P**2 / GEV**2  BETWEEN   0.0   AND  10.                  *
*                       P**2 = 0  <=> REAL PHOTON                 *
*             X         BETWEEN  1.E-4  AND   1.                  *
*                                                                 *
*   HEAVY QUARK THRESHOLDS  Q(H) = M(H)  IN THE BETA FUNCTION :   *
*                   M(C)  =  1.5,  M(B)  =  4.5                   *
*   CORRESPONDING LAMBDA(F) VALUES IN GEV FOR  Q**2 > M(H)**2 :   *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,                                *
*   THE NUMBER OF ACTIVE QUARK FLAVOURS IS  NF = 3  EVERYWHERE    *
*   EXCEPT IN THE BETA FUNCTION, I.E. THE HEAVY QUARKS C,B,...    *
*   ARE NOT PRESENT AS PARTONS IN THE Q2-EVOLUTION.               *
*                                                                 *
*   PLEASE REPORT ANY STRANGE BEHAVIOUR TO :                      *
*          STRAT@HAL1.PHYSIK.UNI-DORTMUND.DE                      *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*...INPUT PARAMETERS :
*
*    X   = MOMENTUM FRACTION
*    Q2  = SCALE Q**2 IN GEV**2
*    P2  = VIRTUALITY OF THE PHOTON IN GEV**2
*
*...OUTPUT (ALWAYS X TIMES THE DISTRIBUTION DIVIDED BY ALPHA_EM) :
*
********************************************************
      subroutine grspar(x,q2,p2,ugam,dgam,sgam,ggam)
      implicit double precision (a-z)
      integer check
c
c     check limits :
c
      check=0
      if(x.lt.0.0001d0) check=1
      if((q2.lt.0.6d0).or.(q2.gt.50000.d0))  check=1
      if(q2.lt.5.d0*p2) check=1
c
c     calculate distributions
c
      if(check.eq.0) then
         call grscalc(x,q2,p2,ugam,dgam,sgam,ggam)
      else
         write(6,*) 'x/q2/p2 - limits exceeded',x,q2,p2
      endif
c
      return
      end
