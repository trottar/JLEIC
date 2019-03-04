*CMZ :  2.08/01 23/06/99  06.42.59  by  Hannes Jung
*CMZ :  2.08/00 06/06/99  16.29.14  by  Hannes Jung
*CMZ :  2.07/03 23/05/99  13.49.18  by  Hannes Jung
*-- Author : Hannes Jung

      FUNCTION PQCDI(X1)
	Implicit None
      REAL X1,PQCDI
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
      REAL XPQ(-6:6)
      REAL SNGL
	Double Precision GEV2NB,ALPH_EM,GF,Y,X,SMIN,CHIMIN,XPMIN,XPMAX
	Double Precision XP,ZQMIN,ZQMAX,XKSI,BETA,BETA2,ZLOG,SH,SCAL
	Double Precision ALPHA_S,ALPHAS
	Double Precision XG,sig0t,sig0s,sigi,sig0,sig,xt,xl,f2,fl,comfac
	Double Precision XQ,wtq,wtg,XOLD
	Integer I,LUCHGE
      DATA GEV2NB/.3893857D+6/
      PQCDI = 0.0
      ALPH_EM = ALPH
      IF(IRUNAEM.EQ.1) ALPH_EM = DBLE(ULALEM(SNGL(Q2)))
      GF = pi*alph_em/sin2w/dsqrt(2.d0)/xmw2
      Y = DBLE(YY)
      X = Q2/Y/SSS
c      write(6,*) ' pqcdi ',x,y,q2,sss
c      write(6,*) 'PQCDI ',x1
      IF(IPRO.NE.14) THEN
         AM(1) = 0.D0
         AM(2) = 0.D0
      ENDIF
      SMIN = DMAX1(4.D0*PT2CUT(IPRO),2.D0*(AM(1)**2+AM(2)**2))
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
c       write(6,*) ' pqcdi ',xpmax,xpmin
c       write(6,*) ' pqcdi ',am(1),am(2),PT2CUT(IPRO),IPRO
c        write(6,*) ' PQCDI: CHIMIN,XPMIN,XPMAX',CHIMIN,XPMIN,XPMAX
         goto 30
      endif
      IF(XPMIN.EQ.0.0D0) THEN
         write(6,*) ' PQCDI XPMIN CHIMIN,X ',CHIMIN,X,IPRO,AM(1)
         write(6,*) ' PQCDI PT2CUT,Y,Q2 ',PT2CUT(IPRO),Y,Q2
      ENDIF
      XP = XPMIN * (XPMAX/XPMIN)**DBLE(X1)
      IF(IPRO.NE.14) THEN
         ZQMIN = 0.5D0 - DSQRT(0.25D0 - PT2CUT(IPRO)*XP/Q2/(1.D0-XP))
         IF(ZQMIN.GT.0.5D0) GOTO 30
         ZQMAX = 1.D0 - ZQMIN
         xksi = 0.D0
         beta =  2.d0*zqmax - 1.d0
      ELSE
         IF((4.d0*am(1)**2*xp/Q2/(1.d0-xp)).GE.1.d0) GOTO 30
         beta2 = 1.d0 - 4.d0*am(1)**2*xp/Q2/(1.d0-xp)
         beta = sqrt(beta2)
         ZQMAX = 1.d0 + beta
         ZQMIN = 1.d0 - beta
         xksi = am(1)**2/q2
c         IF(ZQMIN.GT.0.5) GOTO 123
c         write(6,*) ' pqcdi beta,zqmax,zqmin',beta,zqmax,zqmin
      ENDIF
      zlog = log(zqmax/zqmin)
      sh = Q2/xp*(1.d0 -xp)
	IF(IQ2.EQ.1.AND.AM(1).GT.1.D0) THEN
	   SCAL = 4.D0*AM(1)**2
      ELSEIF(IQ2.EQ.2) THEN
         SCAL = SH
      ELSEIF(IQ2.EQ.4) THEN
         SCAL =4.D0*AM(1)**2 + Q2
      ELSE
         WRITE(6,*) 'PQCD: NO VALID Q2 SCALE. STOP'
	   WRITE(6,*) ' IQ2 = ',IQ2,' mass = ',AM(1)
         STOP
      ENDIF
      SCAL = SCALFA*SCAL
      ALPHA_S=ALPHAS(DSQRT(SCAL))

      IF(IPRO.EQ.13) THEN
         XG = X/XP
         sig0t = ALPHA_s/4.D0/pi/(1.D0-xp)*
     +   (1.D0-xp)*(xp**2+(1.D0-xp)**2)*(2.D0*(zqmin-zqmax)+
     +   log(zqmax**2/zqmin**2))
         sig0s =  ALPHA_S/4.D0/pi*8.D0*xp*(1.D0-xp)*(zqmax-zqmin)
cnew test	
cc	 sig0s = 0.
cend test	
         sigi = 0.D0
         IF(INTER.EQ.2) THEN
            sigi = 0.D0
         ENDIF
         sig0 = 0.5D0*(1.D0+(1.D0-y)**2)*sig0t + (1.D0-y)*sig0s -
     +          ISIGN(1,K(1,2))*y*(1.D0 -0.5D0*y)*sigi
c including phase space factors
         sig = sig0*4.D0*pi*alph_em**2 /y/Q2/Q2
         if(xg.gt.1.) write(6,*) ' pqcdi xg>1 ',xg,ipro,x,xp
         IF(ID.EQ.1) THEN
            CALL PYSTFU(2212,SNGL(XG),SNGL(SCAL),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(XG),SNGL(SCAL),XPQ)
         ENDIF
         wtg = DBLE(XPQ(0))
c         write(6,*) ' 13... Xg,XPQ(0) ',XG,XPQ(0),wtg,X,XP,alpha_s
c          write(6,*) alph,pi
c         write(6,*) y,q2
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
         xt = zlog*(0.5d0*xp - xp**2 *(1.d0-xp) +
     +        2.d0*xksi*xp**2*(1.d0-3.d0*xp) - 4.d0*xksi**2*xp**3)
     +        + beta*(4.d0*xp**2*(1.d0-xp) - 0.5d0*xp
     +        -2.d0*xp**2*xksi*(1.d0-xp))
c note the following line is different in massless case from Peccei/Ruckl
c    +        + beta*(3.d0*xp**2*(1.d0-xp) - 0.5d0*xp ! Peccei/Ruckl
         xl = 2.d0*beta*xp**2*(1.d0-xp) -4.d0*xp**3*xksi*zlog
         f2 = alpha_s * xt /pi
         fl = alpha_s * xl /pi
cnew test	
cc	 fl = 0.
cend test	
c for charged current xf_3 is still missing
         comfac = 4.d0*pi*alph_em**2/y/q2/q2
         sig0 = comfac*(0.5d0*(1.d0+(1.d0-y)**2)*f2 - 0.5d0*y**2*fl)
c including phase space factors
         sig = sig0/xp
         if(xg.gt.1.) write(6,*) ' pqcdi xg>1 ',xg,ipro,x,xp
         IF(ID.EQ.1) THEN
            CALL PYSTFU(2212,SNGL(XG),SNGL(SCAL),XPQ)
         ELSE
            CALL RASTFU(KINT(2,2),SNGL(XG),SNGL(SCAL),XPQ)
         ENDIF
         wtg = dble(XPQ(0))
         sig = sig * wtg
c include charges for charm/bottom
         IF(INTER.LT.2) THEN
            IF(AM(1).LT.2.d0)  sig = sig *4.d0/9.d0
            IF(AM(1).GT.2.d0)  sig = sig *1.d0/9.d0
         ELSEIF(INTER.EQ.2) THEN
            sig = sig * GF**2 * Q2 *Q2/8.D0/pi**2/alph_em**2/(1.d0 + Q2/
     +      XMW2)**2
         ENDIF

      ELSEIF(IPRO.EQ.15) THEN
         sig0t = 2.d0*ALPHA_S/3.d0/pi/(1.d0-xp)*
     +   (xp**2*zlog + zqmin -zqmax +(zqmin**2-zqmax**2)/2.d0 +
     +    zlog+xp*(1.d0-xp)*(zqmax**2-zqmin**2)+
     +    2.d0*(1.d0-xp)*(zqmax-zqmin))
         sig0s = 2.d0*ALPHA_S/3.d0/pi*2*xp*(1.d0-xp)*(zqmax**2-zqmin**2)
     +    /(1.d0-xp)
         sigi = 0.D0
cnew test	
cc	 sig0s = 0.
cend test	
         IF(INTER.EQ.2) THEN
            sigi = 2.d0*ALPHA_S/3.d0/pi/(1.d0-xp)*
     +   (xp**2*zlog + zqmin -zqmax +(zqmin**2-zqmax**2)/2.d0 +
     +    zlog+2.d0*xp*(1.d0-xp)*(zqmax-zqmin)+
     +    2.d0*(1.d0-xp)*(zqmax**2-zqmin**2))
         ENDIF

         XQ=X/XP
         if(xg.gt.1.) write(6,*) ' pqcdi xq>1 ',xq,ipro,x,xp
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
         write(6,*) ' PQCDI process not implemented: IPRO = ',IPRO
      ENDIF
      if(x1.gt.1..or.x1.lt.0.) THEN
         write(6,*) 'PQCDI: x1 = ',x1,IPRO,PT2CUT(IPRO),AM(1)
      endif

      PQCDI = SNGL(SIG * GEV2NB  * XP*LOG(XPMAX/XPMIN))
c      write(6,*) ' pqcdi ',ipro,pqcdi
      RETURN
   30 CONTINUE
      PQCDI = 0.0
      END
