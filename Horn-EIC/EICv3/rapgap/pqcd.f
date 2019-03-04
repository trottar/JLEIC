*CMZ :  2.08/01 23/06/99  07.03.33  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.26.59  by  Hannes Jung
*CMZ :  2.07/03 16/05/99  12.09.03  by  Hannes Jung
*-- Author :    Hannes Jung   28/02/97

      FUNCTION PQCD(X1,Y1)
      Implicit None
      REAL X1,Y1,PQCD
      Double Precision XV,WMAX
      Integer NDIMEN,IMIX
      COMMON /XVAL/ XV(20),NDIMEN
      COMMON /OALPHAS/ WMAX,IMIX
      Integer ID
      double precision xpom
      COMMON /PQPMID/xpom,ID
*KEEP,RGDISDIF.
      INTEGER IDIR,IDISDIF
      COMMON/DISDIF/ IDIR,IDISDIF
*KEEP,RGPARA.
      INTEGER IWEI
      DOUBLE PRECISION ALPHS,PI,ALPH
      COMMON /PARAM/ ALPHS,PI,ALPH,IWEI
      DOUBLE PRECISION SIN2W,XMW2
      COMMON/ELWEAK/SIN2W,XMW2
      DOUBLE PRECISION PT2CUT,THEMA,THEMI,Q2START,W_Q2,OMEG2
      INTEGER IRUNA,IQ2,IRUNAEM
      INTEGER IPRO
      COMMON/RAPA /IPRO,IRUNA,IQ2,IRUNAEM,Q2START,W_Q2,OMEG2
      DOUBLE PRECISION SCALFA
      COMMON/SCALF/ SCALFA
      COMMON/PTCUT/ PT2CUT(100)
      COMMON/ELECT/ THEMA,THEMI
      REAL ULALPS,ULALEM
      EXTERNAL ULALPS,ULALEM
C     SAVE


*KEEP,RGPARA1.
      DOUBLE PRECISION SHAT,YMAX,YMIN,Q2,Q2MAX,Q2MIN
      DOUBLE PRECISION XMAX,XMIN,Q2Q,AM,PCM
      COMMON /PARAT/AM(18),SHAT,YMAX,YMIN,Q2MAX,Q2MIN,XMAX,XMIN
      COMMON /PARAE/Q2,Q2Q,PCM(4,18)
C      SAVE
*KEEP,RGPART.
      DOUBLE PRECISION SSS,CM,DBCMS
      DOUBLE PRECISION PBEAM
      INTEGER KBEAM,KINT
      COMMON /PARTON/ SSS,CM(4),DBCMS(4)
      COMMON /BEAM/PBEAM(2,5),KBEAM(2,5),KINT(2,5)
C      SAVE


*KEEP,RGRAPGKI.
      REAL YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      COMMON/RAPGKI/ YY,XEL,XPR,PT2H,SHH,T2GKI,XFGKI
      REAL ZQGKI,XPGKI,PHITGKI
      COMMON/MEINFO/ZQGKI,XPGKI,PHITGKI
C     SAVE

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

*KEND.
      Double Precision G1,G2,GA,GB,gff,G3,GC
      DIMENSION G1(2),G2(2),GA(2),GB(2)
      DIMENSION gff(4)
      Double Precision H11,H12,H21,H22,H41,H14,H24,H34,H61
      Double Precision ksiw,sg,xg,mf,mfp,lamd,f,g,m2plus,m2minus,m0
      Double Precision Vff,deltal,rho
      Double Precision GEV2NB,ALPH_EM,GF,Y,X,SMIN,CHIMIN,XPMIN,XPMAX
      Double Precision xp,sh,zq,zqmin,zqmax,gam2,PT2,SCAL,t,at
      Double Precision sig0t,sig0s,sigi,sig,sig0,wtg,sum,comfac
      Double Precision xq,wtq,QF2,XOLD
      Double Precision ALPHA_S,ALPHAS
      Integer I,LUCHGE
      REAL XPQ(-6:6)
      REAL SNGL
      DATA GEV2NB/.3893857D+6/
      PQCD = 0.0
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
      GF = pi*alph_em/sin2w/dsqrt(2.d0)/xmw2
      Y = DBLE(YY)
      X = Q2/Y/SSS
      IF(IPRO.NE.14) THEN
         AM(1) = 0.D0
         AM(2) = 0.D0
      ENDIF
      SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
c            write(6,*) ' pqcd  ',x,xpom,id,ipro
      IF(ID.EQ.1) THEN
      ELSE
         XOLD=X
         X=XOLD/XPOM
         IF(X.LT.1.and.X.gt.0.) THEN
c            write(6,*) ' pqcd beta found ',xold,xpom,x
         ELSE
            RETURN
         ENDIF
      ENDIF
cc      CHIMIN=(SMIN+Q2)/(Y*SSS)
      CHIMIN=(1.+SMIN/Q2)*X
      CHIMIN = DMAX1(CHIMIN,X)
      XPMIN = X
      XPMAX = X/CHIMIN
      if(xpmax.lt.xpmin) then
c         write(6,*) 'pqcd: xpmax<xpmin : goto 30'
         goto 30
      endif
      IF(XPMIN.EQ.0.0D0) THEN
         write(6,*) ' PQCD XPMIN CHIMIN,X ',CHIMIN,X,IPRO,AM(1)
         write(6,*) ' PQCD PT2CUT,Y,Q2 ',PT2CUT(IPRO),Y,Q2
      ENDIF
      XP = XPMIN * (XPMAX/XPMIN)**DBLE(X1)
      sh = Q2/xp*(1.d0 -xp)
      IF(IPRO.NE.14) THEN
         ZQMIN = 0.5D0 - DSQRT(0.25D0 - PT2CUT(IPRO)*XP/Q2/(1.D0-XP))
         IF(ZQMIN.GT.0.5D0) THEN
c         write(6,*) 'pqcd: ZQMIN > 0.5 : goto 30'
            GOTO 30
         ENDIF
         ZQMAX = 1.D0 - ZQMIN
      ELSE
         xg = x/xp
         sg = xg*sss
         mf=AM(1)
         mfp=AM(2)
         gam2 = sg*y*(1.d0-xp) + mf**2 - mfp**2
         lamd = (sg*y*(1.d0-xp) - mf**2 - mfp**2)**2 - 4.d0*mf**2*mfp**2
         zqmax = (gam2 + dsqrt(lamd))/sh/2.d0
         zqmin = (gam2 - dsqrt(lamd))/sh/2.d0
      ENDIF
c generate zq and 1-zq pole

      zq = zqmin*(zqmax/zqmin)**dble(y1)
      zq = zqmax/zqmin *(zqmin/zqmax)**(2.d0*dble(y1)) + 1
      zq = 1.d0/zq
      PT2 = sh*zq*(1.d0-zq)
      IF(IQ2.EQ.2) THEN
         SCAL = SH
      ELSEIF(IQ2.EQ.3) THEN
         SCAL = (2.D0*AM(1))**2 + PT2
      ELSEIF(IQ2.EQ.4) THEN
         SCAL = Q2
      ELSEIF(IQ2.EQ.5) THEN
         SCAL = Q2 + PT2 + (2.D0*AM(1))**2
c this is for testing
      ELSEIF(IQ2.EQ.6) THEN
C kt**2 as from Zeppenfeld/Mirkes
         t = -zq*Q2/xp
         at=DSQRT((t+Q2)**2+4.D0*Q2*PT2)
         SCAL = at*(at+t+Q2)/2.D0/Q2
c         write(6,*) 'pqcd : t,PT2,Q2,Q2Q ',t,PT2,Q2,Q2Q
      ELSEIF(IQ2.EQ.7) THEN
c torbjorn suggestion
         SCAL = PT2 + 8.D0*PT2**2/SH

      ELSE
         write(6,*) ' PQCD: IQ2 = ',IQ2
         WRITE(6,*) ' PQCD: NO VALID Q2 SCALE. STOP'
         STOP
      ENDIF
      IF(IQ2.EQ.5) THEN
         SCAL = Q2 + PT2*SCALFA + (2.D0*AM(1))**2
      ELSE
         SCAL = SCALFA*SCAL
      ENDIF
      ALPHA_S=ALPHAS(DSQRT(SCAL))

      IF(IPRO.EQ.13) THEN
         XG = X/XP
         sig0t = ALPHA_s/4.D0/pi/(1.D0-zq)/zq*
     +   (xp**2+(1.D0-xp)**2)*(zq**2+(1.d0-zq)**2)
         sig0s =  ALPHA_S/4.D0/pi*8.D0*xp*(1.D0-xp)
         sigi = 0.D0
         IF(INTER.EQ.2) THEN
            sigi = 0.D0
         ENDIF
         sig0 = 0.5D0*(1.D0+(1.D0-y)**2)*sig0t + (1.D0-y)*sig0s -
     +          ISIGN(1,K(1,2))*y*(1.D0 -0.5D0*y)*sigi
c including phase space factors
         sig = sig0*4.D0*pi*alph_em**2 /y/Q2/Q2
         if(xg.gt.1.) write(6,*) ' pqcd xg>1 ',xg,ipro,x,xp
         IF(ID.EQ.1) THEN
            CALL PYSTFU(2212,SNGL(XG),SNGL(SCAL),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(XG),SNGL(SCAL),XPQ)
         ENDIF
         wtg = DBLE(XPQ(0))
         sig = sig * wtg
c include charges for 3 flavours
         IF(INTER.LT.2) THEN
            sig = sig *2.D0/3.D0
         ELSEIF(INTER.EQ.2) THEN
            sig = sig  * GF**2 * Q2 *Q2/8.D0/pi**2/alph_em**2/(1.d0
     +      + Q2/XMW2)**2
         ENDIF
      ELSEIF(IPRO.EQ.14) THEN
         sig = 0.0d0
         xg = x/xp
         sg = xg*sss
         m2plus = (mf**2 + mfp**2)/y/sg
         m2minus = (mf**2 - mfp**2)/y/sg
         m0 = mf*mfp/y/sg
c         write(6,*) mf,mfp,m0,m2minus,m2plus
c............
         f = (1.D0-xp)*zq*(1.D0-zq) + (zq*(mf**2 - mfp**2) - mf**2)/y/
     +   sg
         g = (1.D0-xp)*(1.D0-zq) + xp*zq + m2minus
         if(f*xp.le.0.) then
c         write(6,*) 'pqcd f,xp,zq,mf**2/y/sg',f,xp,zq,mf**2/y/sg
            return
         endif
         H11 = 4.D0/zq/(1.D0-zq)*(2.D0*zq*(1.D0-zq) + 2.D0*xp*(1.D0-xp)
     +   -1.D0) + 4.D0/zq**2/(1.D0-zq)**2 *(-m2plus**2 + m2plus*
     +   m2minus* (2.D0*zq-1.D0) + m2plus*(2.D0*zq*(zq-1.D0)*(xp-1.D0)
     +   -xp) + m2minus*xp*(2.D0*zq-1.D0))
         H21 = 16.D0*(xp/zq/(1.D0-zq) + (m2plus + (1.D0-zq)*m2minus)/
     +   zq**2/(1D0-zq))
         H41 = 32D0*xp/zq/(1.D0-zq) + 16.D0/zq**2/(1.D0-zq)**2 *
     +   (m2plus + (1.D0-2.D0*zq)*m2minus)
         H61 = - 32.D0*xp/zq/(1.D0-zq) - 32.D0/zq**2/(1.D0-zq) *
     +   (m2plus + (1.D0-zq)*m2minus)
         H12 = - 16D0* m0 *(1.D0-xp)/zq/(1.D0-zq) + 8.D0*m0/zq**2/
     +   (1.D0-zq)**2 * (m2plus + (1.D0-2.D0*zq)*m2minus)
         H22 = - 32.D0*m0/zq/(1.D0-zq)
         H14 = - 8.D0/zq**2/(1.D0-zq)*(-zq**2+zq*xp+m2plus+ (1.D0-zq)*
     +   m2minus)
         H24 = 0.D0
         H34 = 8.D0*(1.D0-2.D0*xp)/zq/(1.D0-zq) - 8.D0*(m2plus+(1.D0-
     +   2.D0*zq)*m2minus)/zq**2/(1.D0-zq)**2
         G1(1) = -2.D0*H11 + f*H41
         G2(1) = 2.D0*xp*f*H41 + H21 + g*(H61 + g * H41)
         GA(1) = DSQRT(xp*f) * (H61 + 2.D0*g*H41)
         GB(1) = f * H41
         G1(2) = - 2.D0* H12
         G2(2) = H22
         GA(2) = 0.d0
         GB(2) = 0.d0
         G3 = H14 + zq*H24 -g*H34
         GC = DSQRT(xp*f)*(H24 - 2.D0*xp*H34)
         IF(INTER.EQ.0) THEN
            IF(AM(1).LT.2.D0) qf2=4.d0/9.d0
            IF(AM(1).GE.2.D0) qf2=1.d0/9.d0
            gff(1) = QF2
            gff(2) = QF2
            gff(4) = 0.D0
            deltal = 0.5d0
         ELSEIF(INTER.EQ.2) THEN
            Vff = 1.d0
            ksiw = 1.D0/(SIN2W*8.D0)*Q2**2/XMW2**2/(1.D0 + Q2/XMW2)**2
            rho = 0.d0
            gff(1) = 4.D0*Vff**2*ksiw**2*(1.D0-ISIGN(1,K(1,2))*rho)
            gff(2) = 0.d0
            gff(4) = 4.D0*Vff**2*ksiw**2*(ISIGN(1,K(1,2)) - rho)
            deltal = 1.d0
         ENDIF
         SUM = gff(1)*(xp*y**2*G1(1) + (1.D0 - y)*G2(1)) +gff(2)
     +   *(xp*y**2*G1(2) + (1.D0 - y)*G2(2)) +gff(4)*(xp*y*(2.D0-y)*G3)
c for charged current xf_3 is still missing
         comfac = ALPHA_S* deltal*2.D0*alph_em**2/y/q2/q2/ 16.D0
         sig = comfac*sum
c including phase space factors
         sig=sig/xp
         if(xg.gt.1.) write(6,*) ' pqcdi xg>1 ',xg,ipro,x,xp
         IF(ID.EQ.1) THEN
            CALL PYSTFU(2212,SNGL(XG),SNGL(SCAL),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(XG),SNGL(SCAL),XPQ)
         ENDIF
         wtg = dble(XPQ(0))
         sig = sig * wtg
      ELSEIF(IPRO.EQ.15) THEN
         sig0t = 2.d0*ALPHA_S/3.d0/pi*(
     +   (zq**2 + xp**2)/(1.d0-xp)/(1.d0-zq) +
     +    2.d0*(xp*zq + 1.d0))
         sig0s = 2.d0*ALPHA_S/3.d0/pi*4.d0*xp*zq
         sigi = 0.D0
         IF(INTER.EQ.2) THEN
            sigi = 2.d0*ALPHA_S/3.d0/pi*(1.d0/(1.d0-xp)/(1.d0-zq)*
     +   (xp**2 + xp**2)+2.d0*(xp+zq))
         ENDIF

         XQ=X/XP
         if(xg.gt.1.) write(6,*) ' pqcd xq>1 ',xq,ipro,x,xp
         IF(ID.EQ.1) THEN
            CALL PYSTFU(2212,SNGL(XQ),SNGL(SCAL),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(XQ),SNGL(SCAL),XPQ)
         ENDIF
         wtq = 0.D0
         IF(INTER.LT.2) THEN
            do 10 I = -NFLQCDC,NFLQCDC
   10       wtq = wtq + DBLE(XPQ(I))*DFLOAT(LUCHGE(I))**2/9.D0
            sig0 = 0.5d0*(1.d0+(1.d0-y)**2)*sig0t + (1.d0-y)*sig0s -
     +      ISIGN(1,K(1,2))*y*(1.D0 -0.5D0*y)*sigi
c including phase space factors
            sig = sig0*4.d0*pi*alph_em**2 /y/Q2/Q2
            sig = sig * wtq
         ELSEIF(INTER.EQ.2) THEN
            sig0 = 0.d0
            DO 20  I=1,NFLAV-1,2
c take only contribution from light quarks
               IF(I.GE.3) GOTO 20
               IF(ISIGN(1,K(1,2)).EQ.1) THEN
                  sig0 = sig0 + 0.5d0*(1.d0+(1.d0-y)**2) *sig0t*
     +            DBLE(XPQ(-I)+XPQ(I+1)) + (1.d0-y)*sig0s*DBLE(XPQ(-
     +            I)+XPQ(I+1)) + ISIGN(1,K(1,2))*y*(1.D0 -0.5D0*y)*
     +            sigi* DBLE(-XPQ(-I)+XPQ(I+1))
               ELSEIF(ISIGN(1,K(1,2)).EQ.-1) THEN
                  sig0 = sig0 + 0.5d0*(1.d0+(1.d0-y)**2) *sig0t*
     +            DBLE(XPQ(I)+XPQ(-I-1)) + (1.d0-y)*sig0s*
     +            DBLE(XPQ(I)+XPQ(-I-1)) - ISIGN(1,K(1,2))*y*(1.D0 -
     +            0.5D0*y)*sigi* DBLE(XPQ(I)-XPQ(-I-1))
               ENDIF

   20       CONTINUE
c including phase space factors
            sig = sig0*4.d0*pi*alph_em**2 /y/Q2/Q2
            sig = sig * GF**2 * Q2 *Q2/8.D0/pi**2/alph_em**2/(1.d0 + Q2/
     +      XMW2)**2
c additional factor 0.5 for handiness of electron
            sig = 0.5D0 * sig
         ENDIF
      ELSE
         write(6,*) ' PQCD process not implemented: IPRO = ',IPRO
      ENDIF
      if(x1.gt.1..or.x1.lt.0.) THEN
         write(6,*) 'PQCD: x1 = ',x1,IPRO,PT2CUT(IPRO),AM(1)
      endif

cccc      PQCD = SNGL(SIG*GEV2NB*XP*LOG(XPMAX/XPMIN)*ZQ*LOG(ZQMAX/ZQMIN))
c use (1-zq) instead of zq because we switched to 1-zq
c      PQCD=SNGL(SIG*GEV2NB*XP*LOG(XPMAX/XPMIN)*(1.-ZQ)*LOG(ZQMAX/ZQMIN))
      PQCD=SNGL(SIG*GEV2NB*XP*LOG(XPMAX/XPMIN)*(1.D0-ZQ)*ZQ*
     +   2.D0*LOG(ZQMAX/ZQMIN))
c       write(6,*) ' pqcd ',ipro,sig,pqcd
      RETURN
   30 CONTINUE
      PQCD = 0.0
c       write(6,*) ' pqcd =0',ipro,sig,pqcd
      END
