! -------------------------------------------------------------------------
! TH (6/2008) - exclusive pi+ generator for EIC
! Generates flat in x, Q2, phi, and calculates cross
! section weight based on Christian Weiss 2008 model.
! Uses standard CERNLIB functions to calculate all 
! kinematics.
! TH (3/21/08) fix problem with xmin>1, should not generate in stot
! TH (4/2008) - add pi0 and klambda production channels
! TH (4/2008) - add missing mass calculation and resolution smearing
! -------------------------------------------------------------------------
! Generate one good event, counting the number of attempts that were made.

	subroutine gen_pipl (nevgen)

	implicit none
#include "mcevnt.inc"
#include "mconsp.inc"
#include "mc_set.inc"
#include "mcEvent.inc"
#include "mcVert.inc"
#include "gen_set.inc"


! TH Define variables - clean up un-needed later

	LOGICAL	ok, lcut

        INTEGER nevgen
	INTEGER charge, gtype, pdgtype
        INTEGER IFAIL
        INTEGER irnd
	INTEGER i,j,k,l,m,n

!	REAL rnd(9)
	REAL rnd(99)
        REAL W12

	REAL myweight,mys ! TH - for testing purposes
        REAL azimphi,sinasym
        REAL ts12,tt,tq2,tdsu,tdsl
	REAL*8 W12d,ttd,eed,ebd,epsd
        REAL r_axis(3), cosphi, sinphi, costheta, sintheta  ! ROTATION TEST
	REAL ppi
	REAL plab3smear,p3res
	REAL plabesmear,peres
	REAL plabsmear
	REAL p4lab4smear,r4lab4smear
	REAL missmass2
	REAL mypt2

! 4-vectors
        REAL v4lab12(4), v4lab2e(4)
        REAL v4lab_e(4), v4lab_eprime(4), v4cm_e(4)
        REAL v4lab1(4), v4lab2(4), v4lab3(4), v4lab4(4)
        REAL v4cm1(4), v4cm2(4), v4cm3(4), v4cm4(4)
        REAL v4cm1e
	REAL v4lab3smear(4),v4lab_eprimesmear(4)
	REAL v4lab4smear(4)
	REAL mmtemp(4),Q2smear,dp,rezest

! 3-vectors
       REAL x_axis(3), y_axis(3), z_axis(3), axis(3,3), v3tmp(3)
       REAL r_hat(3)

! double precision
        REAL*8 m2, m3, m4, s, s12, stot, x, y, nu, Q2, t, tmin
        REAL*8 x_xsec, s12_xsec, W12_xsec
	REAL*8 W2
	REAL*8 cos_thetalab_e, cos_thetacm3, sin_thetacm3
        REAL*8 cos_thetalab, cos_thetacm1

        REAL*8 phi1, phi2, phi3cm
	REAL*8 e1cm, e3cm, p1lab, p1cm, p3cm, omega
        REAL*8 meprime2, dmeprime

	REAL*8 wei,wei3000,regge,sigk
	REAL*8 dsu,dsl

        REAL*8 gam,eps,lminx,lmaxq2,lmint,sigt,sigl,sig

	REAL*8 PION_MASS, PROTON_MASS, NEUTRON_MASS, ELECTRON_MASS

	data ELECTRON_MASS /0.00051099892/
	data PROTON_MASS /0.938271998/
	data NEUTRON_MASS /0.93956533/
	data PION_MASS /0.13957018/

! TH - useful global variables (real)
!	helicity=dble(mcSet_PBValue) 
!	beamcharge

C	mc_excl=.TRUE.

C	if(genSet_reaction.eq.1) then
C  	   write(*,*) 'Using pi+ exclusive reaction',genSet_reaction
C        elseif(genSet_reaction.eq.2) then
C  	   write(*,*) 'Using pi0 exclusive reaction',genSet_reaction
C        elseif(genSet_reaction.eq.3) then
C  	   write(*,*) 'Using K Lambda exclusive reaction',genSet_reaction
C        endif

! Initialize event generation counter

	nevgen = 0

! ============== Generate until one good event is obtained ===============

! Count generated events
! One may want to increment the counter only if the event
! was rejected due to a cut, not bad kinematics.

10	nevgen = nevgen+1


! Simulate an event 
!=================================================================
11	call nran(rnd,20)
        irnd = 1

! TH - define masses using Mandelstam notation

        m2 = PROTON_MASS
        m3 = PION_MASS
	m4 = NEUTRON_MASS
!	p3res=0.004
!	peres=0.004

! TH - Kinematic description of the event vertex, 
! stored in the GEN* variables of /mcevnt/
	if (mcSet_ListEvent) call muprin
	1    ('Generate DIS flat in Q2 and x')

! TH - generate flat in x and Q2, then calculate nu

        genx = mcSet_Xmin+(mcSet_Xmax-mcSet_Xmin)*rnd(irnd)
        irnd = irnd + 1
        x = dble(genx)

! TH - Input for option 1 of mkkine (genq2, gennu, genphi)

        genQ2 = mcSet_Q2Min+(mcSet_Q2Max-mcSet_Q2Min)*rnd(irnd)
        irnd = irnd + 1
        Q2 = dble(genQ2)

        nu = Q2/2./x/m2
        gennu = real(nu)

! TH - check nu cut
       if(gennu.le.mcSet_numin.or.gennu.ge.mcSet_numax) goto 10

! Generate isotropic in azimuthal angle of scattered electron

        phi1 = max(0.,dble(mcSet_PhiMin))
        phi2 = min(twopi,dble(mcSet_PhiMax))
        genphi = real((phi1) + rnd(irnd)*(phi2-phi1))
        irnd = irnd + 1

        call mkkine(1,mcSet_EneBeam,ok)
	if (.not.ok) goto 10


c KINEMATICS ------------------------------------

! TH - DG boost routines somehow give negative energy values - unclear 
! why, for now go with standard CERNLIB

        s12 = (1 - x)/x*Q2 + m2**2 ! gamma*p squared CM energy
        W12 = real(dsqrt(s12))

! TH - target nuleon in the lab frame

        v4lab2(1) = 0.
        v4lab2(2) = 0.
        v4lab2(3) = -mcSet_MomNucleon
        v4lab2(4) = sqrt(mcSet_MomNucleon**2+real(m2)**2)

! TH - beam electron in the lab frame
        v4lab_e(1) = 0.
        v4lab_e(2) = 0.
        v4lab_e(3) = mcSet_EneBeam
        v4lab_e(4) = mcSet_EneBeam

! TH - scattered electron in the lab frame (from mkkine)

	v4lab_eprime(1) = genpx
	v4lab_eprime(2) = genpy
	v4lab_eprime(3) = genpz
	v4lab_eprime(4) = geneprim

! TH - check if mkkine returned a good 4-vector
        meprime2 = geneprim**2 - genpx**2 - genpy**2 - genpz**2
        dmeprime = sqrt(abs(abs(meprime2) - ELECTRON_MASS**2))
        if(dmeprime.ge.(ELECTRON_MASS*2.)) then
!          write(*,*) 'mkkine: electron mass^2 = ', meprime2
	  v4lab_eprime(4) = sqrt(ELECTRON_MASS**2+v4lab_eprime(1)**2
     >                      +v4lab_eprime(2)**2+v4lab_eprime(3)**2)
        endif

	cos_thetalab_e = v4lab_eprime(3)/sqrt(v4lab_eprime(1)**2
     >                 +v4lab_eprime(2)**2+v4lab_eprime(3)**2)


! TH - ep squared CM energy
        call VADD(v4lab2, v4lab_e, v4lab2e, 4)
        stot = v4lab2e(4)**2-v4lab2e(1)**2-v4lab2e(2)**2-v4lab2e(3)**2

! TH - virtual photon in the lab frame
        call VSUB(v4lab_e, v4lab_eprime, v4lab1, 4)


	cos_thetalab = v4lab1(3)/sqrt(v4lab1(1)**2
     >                 +v4lab1(2)**2+v4lab1(3)**2)
	

! TH - Four vector for the CM frame as seen from the lab.
! Motion of the photon + target CM.
        call VADD(v4lab1, v4lab2, v4lab12, 4)

! TH - Boost the photon from lab to CM frame
        call LORENF(W12, v4lab12, v4lab1, v4cm1)

! TH - Boost the beam to the photon/target CM frame
        call LORENF(W12, v4lab12, v4lab_e, v4cm_e)


	cos_thetacm1 = (sqrt(v4cm1(1)**2+v4cm1(2)**2+v4cm1(3)**2)*sqrt(v4cm_e(1)**2+v4cm_e(2)**2+v4cm_e(3)**2))/(v4cm1(1)*v4cm_e(1) + v4cm1(2)*v4cm_e(2) + v4cm1(3)*v4cm_e(3))

! TH - Create a unit vector in the CM frame along the photon momentum 
        call VUNIT(v4cm1,z_axis,3)

! TH - Define the other unit vectors from the cross products
! y_axis = z x e
        call CROSS(z_axis,v4cm_e,y_axis,3)
        call VUNIT(y_axis,y_axis,3)

! TH - x_axis = y x z
        call CROSS(y_axis,z_axis,x_axis,3)
        call VUNIT(x_axis,x_axis,3)


* ==============================================================================

! TH different notation for angles
	phi3cm = dble(rnd(irnd))*PI*2.	! angle of prod plane vs scat plane (phi)
        irnd = irnd + 1

        t = dble(genSet_tmin+(genSet_tmax-genSet_tmin)*rnd(irnd))
        irnd = irnd + 1

! TH - start kinematics tests (using Mandelstam notation for
! photon(1), proton(2), pion(3), and neutron(4)

        omega = (s12+Q2-m2**2)/2./m2
        s = m2**2 + 2*nu*m2 - Q2                  ! not used at the moment
	mys = v4lab12(4)**2-v4lab12(1)**2-v4lab12(2)**2-v4lab12(3)**2 ! check the calculation of s from four vector 

        e1cm = (s12-Q2-m2**2)/2./dble(W12)
        e3cm = (s12+m3**2-m4**2)/2./dble(W12)

        p1lab = dsqrt(Q2+omega**2)

        p1cm = p1lab*m2/dble(W12)

	if (e3cm.gt.m3) then
           p3cm = dsqrt(e3cm**2-m3**2)
	else
	   write(*,*) 'IN EVGEN: E3CM < Mpion',e3cm,m3
	   goto 10
        endif


! TH - cm angles and 4-vector of the meson

	cos_thetacm3 = (t+Q2-m3**2+2.*e1cm*e3cm)/2./p1cm/p3cm

	if(cos_thetacm3**2.ge.1.0) goto 10

        sin_thetacm3 = dsqrt(1.0-cos_thetacm3**2)

! NEW ROTATION 
        cosphi = real(cos(phi3cm))
        sinphi = real(sin(phi3cm))
        costheta = real(cos_thetacm3)
        sintheta = real(sin_thetacm3)
	ppi=real(p3cm)


! TH - comment out for pT2 calculation below
!        call VSCALE(x_axis,cosphi,x_axis,3)
!        call VSCALE(y_axis,sinphi,y_axis,3)
!        call VADD(x_axis,y_axis,r_axis,3)
!        call VSCALE(r_axis,sintheta,r_axis,3)
!        call VSCALE(z_axis,costheta,z_axis,3)
!        call VADD(r_axis,z_axis,v4cm3,3)
!        call VSCALE(v4cm3,ppi,v4cm3,3)
!        v4cm3(4) = real(e3cm)
! TH - end

        call VSCALE(x_axis,cosphi,x_axis,3)
        call VSCALE(y_axis,sinphi,y_axis,3)
        call VADD(x_axis,y_axis,r_hat,3)
        call VSCALE(r_hat,sintheta,r_axis,3)
        call VSCALE(z_axis,costheta,z_axis,3)
        call VADD(r_axis,z_axis,v4cm3,3)
        call VSCALE(v4cm3,ppi,v4cm3,3)
        v4cm3(4) = real(e3cm)
! TH - Calculate pT2
        mypt2=(r_hat(1)*v4cm3(1)+r_hat(2)*v4cm3(2)+r_hat(3)*v4cm3(3))**2 ! full calculation of pt2

! TH Boost pion from CM to lab frame
       call LORENB(W12, v4lab12, v4cm3, v4lab3)
! TH Calculate baryon momentum in lab frame
       call VSUB(v4lab12, v4lab3, v4lab4, 4)


! TH - Apply detector resolutions

! pion
!        call resolution(1,v4lab3,v4lab3smear)
!        call resolution(1,v4lab3,v4lab3smear,dp)
        call resolution(1,v4lab3,v4lab3smear,dp,rezest)


!	write(*,*) 'PION',v4lab3(1),v4lab3(2),v4lab3(3),v4lab3(4),sqrt(v4lab3(4)**2-(v4lab3(1)**2+v4lab3(2)**2+v4lab3(3)**2))
!	write(*,*) 'PION',v4lab3smear(1),v4lab3smear(2),v4lab3smear(3),v4lab3smear(4),sqrt(v4lab3smear(4)**2-(v4lab3smear(1)**2+v4lab3smear(2)**2+v4lab3smear(3)**2))

! electron
!        call resolution(-1,v4lab_eprime,v4lab_eprimesmear)
!        call resolution(-1,v4lab_eprime,v4lab_eprimesmear,dp)
        call resolution(-1,v4lab_eprime,v4lab_eprimesmear,dp,rezest)


!	write(*,*) 'ELECTRONP',v4lab_eprime(1),v4lab_eprime(2),v4lab_eprime(3),v4lab_eprime(4),sqrt(abs(v4lab_eprime(4)**2-(v4lab_eprime(1)**2+v4lab_eprime(2)**2+v4lab_eprime(3)**2)))
!	write(*,*) 'ELECTRONP',v4lab_eprimesmear(1),v4lab_eprimesmear(2),v4lab_eprimesmear(3),v4lab_eprimesmear(4),sqrt(abs(v4lab_eprimesmear(4)**2-(v4lab_eprimesmear(1)**2+v4lab_eprimesmear(2)**2+v4lab_eprimesmear(3)**2)))

! neutron
!        call resolution(0,v4lab4,v4lab4smear)
!        call resolution(0,v4lab4,v4lab4smear,dp)
        call resolution(0,v4lab4,v4lab4smear,dp,rezest)


!	write(*,*) 'NEUTRON',v4lab4(1),v4lab4(2),v4lab4(3),v4lab4(4),sqrt(v4lab4(4)**2-(v4lab4(1)**2+v4lab4(2)**2+v4lab4(3)**2))
!	write(*,*) 'NEUTRON',v4lab4smear(1),v4lab4smear(2),v4lab4smear(3),v4lab4smear(4),sqrt(abs(v4lab4smear(4)**2-(v4lab4smear(1)**2+v4lab4smear(2)**2+v4lab4smear(3)**2)))

!!	write(*,*) 'PX SUM',v4lab3(1)+v4lab_eprime(1)+v4lab4(1)
!!	write(*,*) 'PY SUM',v4lab3(2)+v4lab_eprime(2)+v4lab4(2)
!!	write(*,*) 'PZ SUM',v4lab3(3)+v4lab_eprime(3)+v4lab4(3)

! Calculate missing mass
	call VSUB(v4lab_e,v4lab_eprimesmear,mmtemp,4)
! Calculate Q2 including resolution effects
	Q2smear=mmtemp(4)**2-(mmtemp(1)**2+mmtemp(2)**2+mmtemp(3)**2)
	call VADD(mmtemp,v4lab2,mmtemp,4)
	call VSUB(mmtemp,v4lab3smear,mmtemp,4)

! TH - add this for missing mass resolution checks. Subtracts p_n or p_p
!	call VSUB(mmtemp,v4lab4smear,mmtemp,4)

        missmass2=mmtemp(4)**2-(mmtemp(1)**2+mmtemp(2)**2+mmtemp(3)**2)


!	write(*,*) 'MMTEMP',mmtemp(1),mmtemp(2),mmtemp(3),mmtemp(4),missmass2

!========================================================================
! TH - the cross section models. Note that the inputs are different
! Christian's model needs: W2 (>4), t(<0), Q2(>0) in GeV2 
! sig3000 needs W, Q2, epsilon, -t, sin(thetacm)

! TH - Kinematic checks on x_min and Q2_max, t_min
      CALL LIMX(lminx, stot, Q2, m2)
      if(x.lt.lminx) then
        write(*,*) 'x < x_min',x,lminx
	goto 10
      endif

      lmaxq2 = x*(stot-m2**2)/(1+x*m2**2/(stot-m2**2))
      if(q2.gt.lmaxq2) then
        write(*,*) 'Q2 > Q2_max',Q2,lmaxq2
	goto 10
      endif	

! TH - should already be taken care of above
      CALL FTMIN(stot, Q2, lmint)
      if(abs(t).lt.abs(lmint)) then
         write(*,*) 't < t_min',t,lmint
	 goto 10
      endif	

! TH - Christian's model is only valid for Q2<<W2
      if(Q2.gt.S12) then
!        write(*,*) 'Q2 ~ s',Q2,S12
!	goto 10
      endif	


! TH - need three additional routines by C. Weiss for the function calls
! to FLUX and pip (pip.F etc.)
	W12d=dble(W12)
	eed=dble(v4cm_e(4))
	ttd=abs(t)

! TH - re-calculate W2 for the three different x-ranges
! TH - above x=0.1 use a constant x=0.1 for a given Q2 and t
        if(x.le.0.1) then
           x_xsec=x
           s12_xsec = (1 - x_xsec)/x_xsec*Q2 + m2**2 ! gamma*p squared CM energy
        elseif(x.gt.0.1.and.x.le.0.4) then
   	   x_xsec=0.1
           s12_xsec = (1 - x_xsec)/x_xsec*Q2 + m2**2 ! gamma*p squared CM energy
        else
   	   x_xsec=0.1
           s12_xsec = ((1 - x_xsec)/x_xsec*Q2 + m2**2) ! gamma*p squared CM energy
        endif
        W12_xsec = real(dsqrt(s12_xsec))
	W12d=dble(W12_xsec)


       CALL FLUX(gam, eps, y, stot, x_xsec, Q2, W12d, eed)   ! Photon flux
!       write(*,*) 'TEST WEISS/REGGE FLUX INPUT',stot,x,Q2
!       write(*,*) 'TEST WEISS/REGGE OUTPUT',real(gam),real(eps),y

!          if (s12.gt.4) then
          if (s12_xsec.gt.4) then
!            write(*,*) 'TEST WEISS/REGGE PIP INPUT',s12,t,Q2,x
!             call pip(s12,t,Q2,dsu,dsl) ! Routine by C. Weiss (Regge type model)
! TH - call the function with s12_temp, which is re-calculated value of W2 for
! each of the three x ranges
             call pip(s12_xsec,t,Q2,dsu,dsl) ! Routine by C. Weiss (Regge type model)
!
             sigt = gam*dsu
             sigl = gam*dsl*eps
	     if(x.le.0.4) then
                sig  = sigt + sigl
             else
                sig  = (sigt + sigl)*exp(-x/x_xsec)
             endif
             wei=sig
	     myweight=sig

             regge=dsu+dsl

!      write(*,*) 'TEST WEISS/REGGE OUTPUT',dsu,dsl,sigt,sigl,sig,gam,eps
!      pause
          endif

! TH - Cross section model: Fpi/parm3000 (input W12, )
	  W12d=dble(W12)
	  ttd=abs(t)
!         write(*,*) 'TEST SIG3000 INPUT',W12d,Q2,eps,ttd,sin_thetacm3
!	  pause
	  if(eps.gt.0) then
             call sig3000(W12d,Q2,eps,ttd,sin_thetacm3,wei3000) ! Empirical fit by TH/Fpi12
	  endif
	  wei3000=gam*wei3000 ! include gamma flux factor
      
!========================================================================
! Apply kinematical cuts


! TH - careful with these cuts. They are applied if the logicals are
! set true in the input file. By default: TRUE
! ... to vertex: nu, Q2, x, y, W2
	if (mcSet_cutvertex) then
	   if (mcSet_ListEvent) call muprin ('Checking vertex cuts ...')
	   call mkvertexcuts(lcut)
	   if (lcut) goto 10      ! OK, these events should be counted!
	endif

! ... and to scattered lepton 
	if (mcSet_cutscat) then
	  if (mcSet_ListEvent) call muprin ('Checking scat cuts ...')
          call mktrackcuts(genpprim,genthe,genphi,lcut)
	  if (lcut) goto 10
 	endif

! caculate azimphi
!       azimphi = calcphi(geneprim,genthe,genphi,e_ph,p_ph_x,p_ph_y,p_ph_z)

! If you want to generate sin azimuthal asymmetry, do it here.
        sinasym = 0

! Set the event weight, using the 'weight' variable from /mcevnt/
! (If the event was not generated according to your cross-section,
! the cross-section would be computed here and applied to the weight.
! The result has to come out in microbarns.)

! TH - mcfixedweight calculated in gen_init.F
!	weight = mcfixedweight*wei     
	weight = wei     
!	weight = wei3000

! Store the scattered lepton and any additional tracks in the mcTrack table,
! and store any secondary vertices in the mcVert table.

! TH - start track test

!	write(*,*) '###',beamigtype,beamiltype,beamcharge,beammass,genpx,genpy,genpz
!	pause

! Add the beam lepton to the table mcTrack
      call addmctrack (1,'init',beamigtype,beamiltype,
     $     0,1,1,beamcharge,beammass,real(v4lab_e(1)),
     >    real(v4lab_e(2)),real(v4lab_e(3)))


! ... the target hadron (proton) to the table mcTrack
      charge=1
      gtype=14
      pdgtype=2212
      call addmctrack (2,'init',gtype,pdgtype,
     $     0,1,1,charge,real(m2),real(v4lab2(1)),
     >     real(v4lab2(2)),real(v4lab2(3)))

! ... the virtual photon to the table mcTrack
      charge=0
      gtype=1
      pdgtype=22
      call addmctrack (3,'virt',gtype,pdgtype,
     $     0,1,1,charge,real(v4lab1(4)),real(v4lab1(1)),
     >     real(v4lab1(2)),real(v4lab1(3)))

! ... the scattered beam lepton
      charge=-1
      gtype=3
      pdgtype=11

      call addmctrack (4,'scat',gtype,pdgtype,
     $     0,1,1,charge,beammass,real(v4lab_eprime(1)),
     >     real(v4lab_eprime(2)),real(v4lab_eprime(3)))
C comment out if not including resolution
C     $     0,1,1,charge,beammass,real(v4lab_eprimesmear(1)),
C     >     real(v4lab_eprimesmear(2)),real(v4lab_eprimesmear(3)))


! ... the scattered pion
      charge=1
      gtype=8
      pdgtype=211
      call addmctrack (5,'elas',gtype,pdgtype,
     $     0,1,0,charge,real(m3),real(v4lab3(1)),
     >     real(v4lab3(2)),real(v4lab3(3)))
C comment out if not including resolution
C     $        0,1,0,charge,real(m3),real(v4lab3smear(1)),
C     >        real(v4lab3smear(2)),real(v4lab3smear(3)))

! ... the recoiling neutron
      charge=0
      gtype=13
      pdgtype=2112
      call addmctrack (6,'elas',gtype,pdgtype,
     $     0,1,0,charge,real(m4),real(v4lab4(1)),
     >     real(v4lab4(2)),real(v4lab4(3)))
C comment out if not including resolution
C     $        0,1,0,charge,real(m4),real(v4lab4smear(1)),
C     >        real(v4lab4smear(2)),real(v4lab4smear(3)))
 

! TH - end track test

! =============== End generation section, we're happy ===============

! Fill the mcEvent table

! ... these entries have to be set by hand, depending on your generator
! ........ a code indicating the process which produced this event
	mcEvent_process = 0
! ........ the struck nucleon 
! ........ (use LEPTO convention: 1=proton, 2=neutron)
	mcEvent_nucleon = 0
! ........ the flavour of the struck nucleon
! ........ (use Lund/LEPTO convention: 1=d, 2=u, 3=s, 4=c, -1=dbar, ...)
	mcEvent_qflavour = 0
! ........ index of struck quark in mcTrack table
	mcEvent_qtrack = 0
! ........ a code indicating what happened to the target remnant
	mcEvent_remnant = 0

! ... the remaining entries can be filled automatically, 
! ... just provide a 4-character string for mcEvent_cType
	call fillmcevent('EXCL')

! Fill the mcUser table with any generator-specific event information
! you would like to send to the output file.
! The name of the entry should be at most 16 characters long,
! and the value should be a real number

! TH - to read this out edit file exclusive.in
!	call addmcuser(1,'MyVariable',1.)
	call addmcuser(1,'cos_thetacm3',real(cos_thetacm3))
	call addmcuser(2,'t',real(t))
	call addmcuser(3,'missmass2',real(missmass2))
	call addmcuser(4,'e_n_lab',real(v4lab4(4)))
  	call addmcuser(5,'px_pi_lab',real(v4lab3(1)))
	call addmcuser(6,'py_pi_lab',real(v4lab3(2)))
	call addmcuser(7,'pz_pi_lab',real(v4lab3(3)))
	call addmcuser(8,'px_n_lab',real(v4lab4(1)))
	call addmcuser(9,'py_n_lab',real(v4lab4(2)))
	call addmcuser(10,'pz_n_lab',real(v4lab4(3)))
        call addmcuser(11,'dp',real(dp))
	call addmcuser(12,'px_pi_smear',real(v4lab3smear(1)))
	call addmcuser(13,'py_pi_smear',real(v4lab3smear(2)))
	call addmcuser(14,'pz_pi_smear',real(v4lab3smear(3)))
!	call addmcuser(12,'px_pi_smear',real(v4lab_eprimesmear(1)))
!	call addmcuser(13,'py_pi_smear',real(v4lab_eprimesmear(2)))
!	call addmcuser(14,'pz_pi_smear',real(v4lab_eprimesmear(3)))
!	call addmcuser(15,'q2smear',real(Q2smear))
	call addmcuser(15,'epsilon',real(eps))
        call addmcuser(16,'stot',real(stot))
!        call addmcuser(16,'rezest',real(rezest))
!	call MUPRIN('good event.')

! Increase the value here - more random numbers in the array above
        if(irnd.gt.99) call muprin('irnd > 99')

!	write(*,*) 'MADE IT'

	return
	end

! =================================================================


