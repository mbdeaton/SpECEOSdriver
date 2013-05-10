program driver

  use eosmodule
  implicit none

  interface
	FUNCTION rtnewt_findYe(funcd,x1,x2,xacc,xrho,xtemp,tableymin)
	!USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xtemp,tableymin
	REAL*8 :: rtnewt_findYe
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv,xrho,xtemp,tableymin)
		!USE nrtype
		IMPLICIT NONE
		REAL*8, INTENT(IN)  :: x,xrho,xtemp,tableymin
		REAL*8, INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	END FUNCTION rtnewt_findYe 
	FUNCTION rtnewt_findtheta(funcd,x1,x2,xacc,xrho,xye,xent)
	!USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xye,xent
	REAL*8 :: rtnewt_findtheta
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv,xrho,xye,xentGoal)
		!USE nrtype
		IMPLICIT NONE
		REAL*8, INTENT(IN)  :: x,xrho,xye,xentGoal
		REAL*8, INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	END FUNCTION rtnewt_findtheta
	SUBROUTINE mu_mismatch(x,fval,fderiv,xrho,xtemp,tableymin)
	!USE nrtype
	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: x,xrho,xtemp,tableymin
	REAL*8, INTENT(OUT) :: fval,fderiv
	END SUBROUTINE mu_mismatch
	SUBROUTINE ent_mismatch(x,fval,fderiv,xrho,xye,xentGoal)
	!USE nrtype
	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: x,xye,xrho,xentGoal
	REAL*8, INTENT(OUT) :: fval,fderiv
	END SUBROUTINE ent_mismatch
  end interface

  real*8 xrho,xye,xtemp,tableymin
  real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
  real*8 xdpderho,xdpdrhoe
  integer keytemp,keyerr

  ! for full eos call:
  real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  real*8 xxa,xxh,xxn,xxp

  integer i,j,k
  integer nr, nt, ny
  real*8 ltmax, ltmin, lymin, lymax, lrmin, lrmax
  real*8 lrmin_for_user_chooses_low_rho_bound
  real*8 tmax, tmin, ymin, ymax, rmin, rmax, tminout, tmaxout
  real*8 xunit,tunit,munit,rhounit,punit,eomunit,c2unit
  real*8 GammaP, EoverRho
  real*8 RhoScale, HeatCapacity, GammaTh, kappa, gtfac
  real*8 BARYMASS, RTACC
  real*8 r_tmin, r_tmax
  real*8 theta, xrhothet, xentthet
  real*8 pgam,delta_t,delta_y,xtemp_p,xtemp_pp,xye_p
  real*8 xprs_tp,xent_tp,xenr_tp,xprs_tpp,xent_tpp,xenr_tpp
  real*8 xprs_yp,xent_yp,dPdT,dsdT,dPdY,dsdY,dPdYs,dPds
  real*8 hmin,cs
  real*8 rho_at_hmin,temp_at_hmin,ye_at_hmin
  integer FULL,BETA,COLD,MICRO,DERIVS,MUX,COLDYE,BETAYE,eostype
  integer HSHEN_2011,LS_220,GSHEN_NL3,GSHEN_FSU21,HEMPEL_SFHO,HEMPEL_SFHX,HEMPEL_DD2,eos
  logical USER_CHOOSES_BOUNDS
  logical USER_CHOOSES_LOW_RHO_BOUND

  RTACC = 1.e-5
  keytemp = 1
  keyerr  = 2
  RhoScale = 5.e-5   ! this currently isn't used

  ! Conversion from cgs to geometric (+M_sol=1)
  !   Note: in this code, rho is in g/cc, t in MeV. In SpEC tabulated EOSs we use geo
  !   lrho_cgs = 17.79 + lrho_geo
  xunit = 1.4766d5
  tunit = 4.9255d-6
  munit = 1.989d33
  BARYMASS = 1.659d-24
  rhounit = munit/xunit**3
  punit = munit/(xunit*tunit**2)
  c2unit = (xunit/tunit)**2
  eomunit = c2unit

  ! Define types of tables:
  ! Tables of full thermodynamic potentials
  FULL = 1     ! internal energy in 3D table
  BETA = 2     ! internal energy in 2D table, with ye=beta_equil
  COLD = 3     ! internal energy in 1D table, with t=tmin, ye=beta_equil
  ! Tables of ye
  BETAYE = 4   ! Ye in 2D table, with ye=beta_equil
  COLDYE = 5   ! Ye in 1D table, with t=tmin, ye=beta_equil
  ! Tables of other thermodynamic potentials
  MICRO = 6    ! microphysical potentials in 3D table
  MUX = 7      ! chemical potentials in 3D table
  DERIVS = 8   ! thermodynamic derivs in 3D table (for testing)

  ! Define versions of tables:
  HSHEN_2011 = 1   ! HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5
  LS_220 = 2       ! LS220_234r_136t_50y_analmu_20091212_SVNr26.h5
  GSHEN_NL3 = 3    ! GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5
  GSHEN_FSU21 = 4  ! GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5
  HEMPEL_SFHO = 5  ! Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5
  HEMPEL_SFHX = 6  ! Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5
  HEMPEL_DD2 = 7   ! Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5

  ! ***** User-Chosen Parameters ************************************************************
  eostype = BETAYE
  eos = HEMPEL_DD2

  USER_CHOOSES_LOW_RHO_BOUND = .true. ! true if user is to override table's low bound in r
  ! Choose low density bound different than table's intrinsic bound in r.
  !   This option is often used rather than USER_CHOOSES_BOUNDS because we want to honor
  !   the table's intrinsic bounds in all vars, except the minimum rho, where we do not need
  !   these lowest densities because we apply atmosphere treatment near rho=1e-10.
  !   if USER_CHOOSES_LOW_RHO_BOUND then user must supply low density bound here
  !   if .not.USER_CHOOSES_LOW_RHO_BOUND then this is never used
  lrmin_for_user_chooses_low_rho_bound = 8d0

  USER_CHOOSES_BOUNDS = .false. ! true if user is to override the table's bounds in r,t,y
  ! Choose bounds different than table's intrinsic bounds in r,t,y.
  !   This option is rarely used because the table's intrinsic bounds are representative
  !   of the range of expected r,t,y achieved in evolutions.
  !   if USER_CHOOSES_BOUNDS then user must supply all of the bounds here
  !   if .not.USER_CHOOSES_BOUNDS then these are all overwritten
  lrmin = 8d0
  lrmax = 15.3d0
  ltmin = -2d0
  ltmax = 1.875d0
  ymin = 0.05d0
  ymax = 0.559d0

  if (USER_CHOOSES_BOUNDS.and.USER_CHOOSES_LOW_RHO_BOUND) then
    write(*,"(A,I1,A)") 'Error, user has specified a lower rho bound twice.'
  end if

  ! choose output resolution (nt and/or ny overwritten for some eostypes)
  !nr = 2000
  nr = 250
  nt = 120
  ny = 100

  ! ***** reset table limits and resolution *****************************************************
  if (eos.eq.HSHEN_2011) then
    call readtable("HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5")
    tableymin = 0.01d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 3d0
      lrmax = 16d0
      ltmin = -2.0d0
      ltmax = 2.5d0
      ymin = tableymin
      ymax = 0.6499d0 ! the largest ye in the *h5 table is 0.65, but that makes rootfinding fail here
    end if
  else if (eos.eq.LS_220) then
    call readtable("LS220_234r_136t_50y_analmu_20091212_SVNr26.h5")
    tableymin = 0.035d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 3d0
      lrmax = 16d0
      ltmin = -2.0d0
      ltmax = 2.4d0
      ymin = tableymin
      ymax = 0.53d0
    end if
  else if (eos.eq.GSHEN_NL3) then
    call readtable("GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5")
    tableymin = 0.05d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 3d0
      lrmax = 15.39d0
      ltmin = -2d0
      ltmax = 2.5d0
      ymin = tableymin
      ymax = 0.5599d0 ! the largest ye in the *h5 table is 0.56, but that makes rootfinding fail here.
    end if
  else if (eos.eq.GSHEN_FSU21) then
    call readtable("GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5")
    tableymin = 0.05d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 3d0
      lrmax = 15.39d0
      ltmin = -2d0
      ltmax = 2.5d0
      ymin = tableymin
      ymax = 0.5599d0 ! the largest ye in the *h5 table is 0.56, but that makes rootfinding fail here.
    end if
  else if (eos.eq.HEMPEL_SFHO) then
    call readtable("Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5")
    tableymin = 0.01d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 2.22025d0
      lrmax = 15.5002d0
      ltmin = -2d0
      ltmax = 2.2d0
      ymin = tableymin
      ymax = 0.599d0 ! the largest ye in the *h5 table is 0.6, but that makes rootfinding fail here.
    end if
  else if (eos.eq.HEMPEL_SFHX) then
    call readtable("Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5")
    tableymin = 0.01d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 2.22025d0
      lrmax = 16.2202d0
      ltmin = -2d0
      ltmax = 2.2d0
      ymin = tableymin
      ymax = 0.599d0 ! the largest ye in the *h5 table is 0.6, but that makes rootfinding fail here.
    end if
  else if (eos.eq.HEMPEL_DD2) then
    call readtable("Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5")
    tableymin = 0.01d0
    if (.not.USER_CHOOSES_BOUNDS) then
      lrmin = 2.22025d0
      lrmax = 16.2202d0
      ltmin = -2d0
      ltmax = 2.2d0
      ymin = tableymin
      ymax = 0.599d0 ! the largest ye in the *h5 table is 0.6, but that makes rootfinding fail here.
    end if
  else
    write(*,"(A,I1,A)") 'Error, eos ', eos, ' is not recognized.'
  end if

  ! reset lrmin to user-chosen bound if this option is used
  if (USER_CHOOSES_LOW_RHO_BOUND) then
     lrmin = lrmin_for_user_chooses_low_rho_bound
  end if

  rmin = 10**lrmin
  rmax = 10**lrmax
  tmin = 10**ltmin
  tmax = 10**ltmax
  lymin = log10(ymin)
  lymax = log10(ymax)

  ! reset output resolutions for lower-dimensional tables
  if((eostype.eq.COLD).or.(eostype.eq.COLDYE)) then
     nt = 1
     ny = 1
  end if
  if((eostype.eq.BETA).or.(eostype.eq.BETAYE)) then
     ny = 1
  end if

  ! ***** Print Header ************************************************************************
  if(eostype.eq.MICRO) then
    write(*,"(A)") "EoSType = Microphysics"
  else if(eostype.eq.MUX) then
    write(*,"(A)") "EosType = ChemicalPotentials"
  else if(eostype.eq.DERIVS) then
    write(*,"(A)") "EoSType = Derivs"
  else if(eostype.eq.BETA) then
    write(*,"(A)") "# 2D beta-equilibrium table"
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.FULL) then
    write(*,"(A)") "# Full 3D table"
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.COLD) then
    write(*,"(A,E15.6)") "# Cold beta-equilibrium table, T = ", tmin
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.COLDYE) then
    write(*,"(A,E15.6)") "# T = ", tmin
    write(*,"(A)") "EoSType = ColdBetaYe"
  else if(eostype.eq.BETAYE) then
    write(*,"(A)") "# 2D beta-equilibrium table"
    write(*,"(A)") "EoSType = HotBetaYe"
  else
    write(*,"(A,I1,A)") 'Error, eostype ', eostype, ' is not recognized.'    
  end if

  if(eostype.eq.COLDYE) then
    write(*,"(A,I6)") "Nrho =",nr
  else if(eostype.eq.BETAYE) then
    write(*,"(A,I6,A,A,I6)") "Nrho =",nr,"     ","NT =",nt
  else  
    write(*,"(A,I6,A,A,I6,A,A,I6)") "Nrho =",nr,"     ","Nye =",ny,"     ","NT =",nt
  end if

  write(*,"(A, E15.6, A, A, E15.6)") "RhoMin =",rmin/rhounit,"      ",&
        "RhoMax =",rmax/rhounit

  if((eostype.eq.MICRO).or.(eostype.eq.FULL).or.(eostype.eq.DERIVS) &
    .or.(eostype.eq.MUX)) then
    write(*,"(A, E15.6, A, A, E15.6)") "YeMin =",ymin,"      ",&
      "YeMax =",ymax
  end if

  if((eostype.eq.MICRO).or.(eostype.eq.FULL).or.(eostype.eq.BETA) &
    .or.(eostype.eq.DERIVS).or.(eostype.eq.MUX).or.(eostype.eq.BETAYE)) then
    write(*,"(A, E15.6, A, A, E15.6)") "Tmin =",tmin,"      ",&
      "Tmax =",tmax
  end if

  kappa = 100.d0    ! for our output of gamma s.t. pressure=kappa*rho^gamma
  HeatCapacity = 1.6d-3
  GammaTh = 5.d0/3.d0
  gtfac = (GammaTh-1.d0)/GammaTh
  if((eostype.ne.MICRO).and.(eostype.ne.DERIVS).and.(eostype.ne.MUX) &
	.and.(eostype.ne.COLDYE).and.(eostype.ne.BETAYE)) then
    write(*,"(A, E15.6)") "HeatCapacity =",HeatCapacity,&
      "GammaTh =",GammaTh,&
      "Kappa =",kappa
  end if
  write(*,"(A)") "RhoSpacing = Log"
  if((eostype.eq.MICRO).or.(eostype.eq.FULL).or.(eostype.eq.DERIVS) &
	.or.(eostype.eq.MUX)) then
    write(*,"(A)") "YeSpacing = Linear"
  end if
  if((eostype.eq.MICRO).or.(eostype.eq.FULL).or.(eostype.eq.BETA) &
	.or.(eostype.eq.DERIVS).or.(eostype.eq.MUX).or.(eostype.eq.BETAYE)) then
    write(*,"(A)") "TSpacing = Log"
  end if

  ! ***** Print Table ************************************************************************
  ! Here we print the input table used in a SpEC Hydro evolution.
  !   To be read in using the EquationOfState::Tabulated class.
  !   For the standard tables -- FULL (3D), BETA (2D), COLD(1D) -- as of 3.5.12,
  !   the output colums are:
  ! [1] epsilon, specific internal energy
  ! [2] gamma,   so that P=kappa*rho^gamma, where kappa is the constant defined above
  ! [3] cs,      relativistic adiabatic sound speed

  ! In the current implementation hmin vars are not used; they were used once to find hmin in the table.  
  hmin = 1.d0
  rho_at_hmin = rmin
  temp_at_hmin = tmin
  ye_at_hmin = ymin
  do i=0, nt-1
    if((eostype.eq.COLD).or.(eostype.eq.COLDYE)) then
      xtemp = tmin
    else
      !xtemp = tmin + i*(tmax-tmin)/nt
      xtemp = 10**(ltmin + i*(ltmax-ltmin)/(nt-1))
    end if
    do j=0, ny-1
      !xye = 10**(lymin + j*(lymax-lymin)/(ny-1))
      xye = ymin + j*(ymax-ymin)/(ny-1)
      do k=0, nr-1
        xrho = 10**(lrmin + k*(lrmax-lrmin)/(nr-1))
	!xrho  = 1.d-3*rhounit
	!xtemp = 8.19408d0
        !xtemp = 1.6d0
	!xye   = 0.15d0
	! set xye to satisfy beta_equil
        if((eostype.eq.BETA).or.(eostype.eq.COLD).or.(eostype.eq.COLDYE) &
	.or.(eostype.eq.BETAYE)) then
          xye = rtnewt_findYe(mu_mismatch,ymin,ymax,RTACC,xrho,xtemp,tableymin)
        end if

	! Extrapolate to T=0
	!   WARNING: careful here, I've moved pieces of code around without thinking about
	!   dependencies (Brett)
	!if((eostype.eq.COLDYE).or(eostype.eq.COLD)) then
	!  xtemp_p = 2.d0*xtemp
	!  call nuc_eos_full(xrho,xtemp_p,xye,xenr_tp,xprs_tp,xent_tp,xcs2, &
	!    xdedt,xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
	!    xmu_e,xmu_n,xmu_p,xmuhat,keytemp,keyerr)
	!  xtemp_pp = 3.d0*xtemp
	!  call nuc_eos_full(xrho,xtemp_pp,xye,xenr_tpp,xprs_tpp,xent_tpp, &
	!    xcs2, xdedt,xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
	!    xmu_e,xmu_n,xmu_p,xmuhat,keytemp,keyerr)
	!end if
	!xenr_tp = xenr_tp/eomunit
	!xprs_tp = xprs_tp/punit
	!xenr_tpp = xenr_tpp/eomunit
	!xprs_tpp = xprs_tpp/punit
	!xenr = 3.d0*xenr - 3.d0*xenr_tp + xenr_tpp
	!xprs = 3.d0*xprs - 3.d0*xprs_tp + xprs_tpp 

        call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)

        if(eostype.eq.MICRO) then
	  xrhothet = 3.e12
          theta = rtnewt_findtheta(ent_mismatch,tmin,tmax,RTACC,&
	    xrhothet,xye,xent)
        end if

	if(eostype.eq.DERIVS) then	
	  pgam = xcs2*xrho
	  delta_t = xtemp*0.01
	  delta_y = xye*0.01
	  xtemp_p = xtemp + delta_t
	  xye_p   = xye + delta_y
          call nuc_eos_full(xrho,xtemp_p,xye,xenr,xprs_tp,xent_tp,xcs2,xdedt,&
            xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
            xmuhat,keytemp,keyerr)
          dPdT = (xprs_tp - xprs)/delta_t
	  dsdT = (xent_tp - xent)/delta_t
          call nuc_eos_full(xrho,xtemp,xye_p,xenr,xprs_yp,xent_yp,xcs2,xdedt,&
            xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
            xmuhat,keytemp,keyerr)
          dPdY = (xprs_yp - xprs)/delta_y
	  dsdY = (xent_yp - xent)/delta_y

	  dPdYs = dPdY + dPdT*(dsdY/dsdT)
	  dPds  = dPdT/dsdT

	  dPdY = dPdY/punit
	  dPds = dPds/punit
	  pgam = pgam/(c2unit*rhounit)
	end if

        xrho = xrho/rhounit
        xenr = xenr/eomunit
        xcs2 = xcs2/c2unit
        xprs = xprs/punit
	xent = xent*9.56565348d17/eomunit

        GammaP = log(xprs/kappa)/log(xrho)
        EoverRho = xenr

        if(xxn<1.e-20) xxn = 0.d0
        if(xxp<1.e-20) xxp = 0.d0
        if(xxa<1.e-20) xxa = 0.d0
        if(xxh<1.e-20) xxh = 0.d0

	if(1.0+xenr+xprs/xrho.lt.hmin) then
	  hmin = 1.0+xenr+xprs/xrho
	  rho_at_hmin = xrho
	  temp_at_hmin = xtemp
	  ye_at_hmin = xye
	end if

        if(eostype.eq.MICRO) then
          write(*,"(1P10E15.6)") xent,xxn,xxp,xxa,xxh,xabar,xzbar,theta	
	else if(eostype.eq.MUX) then
	  write(*,"(1P10E15.6)") xmu_e,xmu_p-xmu_n,xxn,xxp,xxh,xabar,xzbar
	else if(eostype.eq.DERIVS) then
	  write(*,"(1P10E15.6)") pgam,dPds,dPdY,xent
	else if(eostype.eq.COLDYE) then
	  write(*,"(1P10E15.6)") xrho, xye
	else if(eostype.eq.BETAYE) then
	  write(*,"(1P10E15.6)") xye
        else if ((eostype.eq.FULL).or.(eostype.eq.BETA).or.(eostype.eq.COLD)) then
	  cs = 0.d0
	  if(xcs2.gt.0.d0) cs = sqrt(xcs2/(1.0+xenr+xprs/xrho))
          write(*,"(1P10E15.8)") EoverRho,GammaP,cs
	else
          write(*,"(A,I1,A)") 'Error, eostype ', eostype, ' is not recognized.'    
        end if
      end do ! END for r
    end do ! END for y
  end do ! END for t

  ! MinimumEnthalpy can be temporarily written to the end of the table for examination
  !write(*,"(A,E15.6)") "MinimumEnthalpy = ",hmin
  !write(*,"(A,E15.6)") "rho_at_hmin = ",rho_at_hmin
  !write(*,"(A,E15.6)") "temp_at_hmin = ",temp_at_hmin
  !write(*,"(A,E15.6)") "ye_at_hmin = ",ye_at_hmin

end program driver

! ****************************************************************************************
SUBROUTINE mu_mismatch(xye,fval,fderiv,xrho,xtemp,tableymin)
  !USE nrtype
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xye,xrho,xtemp,tableymin
  REAL*8, INTENT(OUT) :: fval,fderiv
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: ye_p, ye_m, EPS, f_p, f_m, deriv_crit, ymin, lymin
  EPS = 1.d-7
  keytemp = 1
  keyerr  = 0

  !write(*,*) "1finding chemical potentials ",xye,xrho,xtemp
  call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  fval = xmu_n - xmu_p - xmu_e
  !write(*,*) "fval=", fval," n=",xmu_n,", p=", xmu_p,", e=",xmu_e

  ! set lower values for symmetric deriv
  ye_m = xye - EPS
  if(ye_m.lt.tableymin) ye_m = tableymin
  call nuc_eos_full(xrho,xtemp,ye_m,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  f_m = xmu_n - xmu_p - xmu_e

  ! set upper values for symmetric deriv
  ye_p = ye_m + 2.d0*EPS
  call nuc_eos_full(xrho,xtemp,ye_p,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  f_p = xmu_n - xmu_p - xmu_e

  fderiv = (f_p - f_m)/(2.d0*EPS)

  ! keep next guess within ye-bounds
  deriv_crit = 1.01d0*fval/(xye-tableymin)
  if(abs(fderiv).lt.abs(deriv_crit)) then
    fderiv = deriv_crit
  end if
END SUBROUTINE mu_mismatch

FUNCTION rtnewt_findYe(funcd,x1,x2,xacc,xrho,xtemp,tableymin)
  !USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xtemp,tableymin
  REAL*8 :: rtnewt_findYe
  INTERFACE
	SUBROUTINE funcd(x,fval,fderiv,xrho,xtemp,tableymin)
	  !USE nrtype
	  IMPLICIT NONE
	  REAL*8, INTENT(IN) :: x,xrho,xtemp,tableymin
	  REAL*8, INTENT(OUT) :: fval,fderiv
	END SUBROUTINE funcd
  END INTERFACE
  INTEGER, PARAMETER :: MAXIT=40
  INTEGER :: j
  REAL*8 :: df,dx,f,fl,fh,dfl,dfh, xl,xh
  !rtnewt_findYe=0.5d0*(x1+x2)
  !call funcd(rtnewt_findYe,f,df,xrho,xtemp,tableymin)
  if(x1<x2) then
    xl = x1
    xh = x2
  else
    xl = x2
    xh = x1
  end if
  call funcd(xl,fl,dfl,xrho,xtemp,tableymin)
  call funcd(xh,fh,dfh,xrho,xtemp,tableymin)

  !write(*,*) "fl=",fl,"fh=",fh

  if(abs(fl)<xacc) then
    rtnewt_findYe = xl
    RETURN
  end if
  if(abs(fh)<xacc) then
    rtnewt_findYe = xh
    RETURN
  end if

	if(fl*fh>0.d0) then
	  rtnewt_findYe = xl
	  RETURN
	else 
	  rtnewt_findYe=0.5d0*(xl+xh)
	  do j=1,MAXIT
	    call funcd(rtnewt_findYe,f,df,xrho,xtemp,tableymin)
	    if(abs(f)<xacc) RETURN
	    if(f*fl<0.d0) then
		xh = rtnewt_findYe
		fh = f
	    else
		xl = rtnewt_findYe
		fl = f
	    end if
	    rtnewt_findYe = 0.5d0*(xl+xh)
	  end do
	end if
!	do j=1,MAXIT
!		call funcd(rtnewt_findYe,f,df,xrho,xtemp,tableymin)
!		dx=f/df
!		rtnewt_findYe=rtnewt_findYe-dx
!		if ((x1-rtnewt_findYe)*(rtnewt_findYe-x2) < 0.0)&
!                     write(*,*) "rtnewt_findYe: values jumped out of brackets"
!			!call nrerror('rtnewt_findYe: values jumped out of brackets')
!		if (abs(dx) < xacc) RETURN
!	end do
	!call nrerror('rtnewt_findYe exceeded maximum iterations')
        write(*,*) "rtnewt_findYe exceeded maximum iterations"
END FUNCTION rtnewt_findYe

SUBROUTINE ent_mismatch(xtemp,fval,fderiv,xrho,xye,xentGoal)
  !USE nrtype
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xrho,xye,xtemp,xentGoal
  REAL*8, INTENT(OUT) :: fval,fderiv
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: temp_p, temp_m, EPS, f_p, f_m, deriv_crit, tmin
  EPS = 1.d-7
  keytemp = 1
  keyerr  = 0


!  xrho = 3.d12

  call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  fval = xent - xentGoal

!  temp_p = xye + EPS
!  call nuc_eos_full(xrho,temp_p,xye,xenr,xprs,xent,xcs2,xdedt,&
!       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
!       xmuhat,keytemp,keyerr)
!  f_p = xent - xentGoal

!  temp_m = xtemp - EPS
!  call nuc_eos_full(xrho,temp_p,xye,xenr,xprs,xent,xcs2,xdedt,&
!       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
!       xmuhat,keytemp,keyerr)
!  f_m = xent - xentGoal

!  fderiv = (f_p - f_m)/(2.d0*EPS)
  tmin = 0.1d0
  deriv_crit = 1.01d0*fval/(xtemp-tmin)
  if(abs(fderiv).lt.abs(deriv_crit)) then
    fderiv = deriv_crit
  end if
END SUBROUTINE ent_mismatch

	FUNCTION rtnewt_findtheta(funcd,x1,x2,xacc,xrho,xye,xent)
	!USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xye,xent
	REAL*8 :: rtnewt_findtheta
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv,xrho,xye,xentGoal)
		!USE nrtype
		IMPLICIT NONE
		REAL*8, INTENT(IN) :: x,xye,xrho,xentGoal
		REAL*8, INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	INTEGER, PARAMETER :: MAXIT=40
	INTEGER :: j
	REAL*8 :: df,dx,f,fl,fh,dfl,dfh, xl,xh
!	rtnewt_findtheta=0.5d0*(x1+x2)
!	call funcd(rtnewt_findtheta,f,df,xrho,xye,xent)
	if(x1<x2) then
	  xl = x1
	  xh = x2
	else
	  xl = x2
	  xh = x1
	end if
	call funcd(xl,fl,dfl,xrho,xye,xent)
	call funcd(xh,fh,dfh,xrho,xye,xent)
	if(abs(fl)<xacc) then
	  rtnewt_findtheta = xl
	  RETURN
	end if
	if(abs(fh)<xacc) then
	  rtnewt_findtheta = xh
	  RETURN
	end if

	if(fl*fh>0.d0) then
	  rtnewt_findtheta = xl
	  RETURN
	else 
	  rtnewt_findtheta=0.5d0*(xl+xh)
	  do j=1,MAXIT
	    call funcd(rtnewt_findtheta,f,df,xrho,xye,xent)
	    if(abs(f)<xacc) RETURN
	    if(f*fl<0.d0) then
		xh = rtnewt_findtheta
		fh = f
	    else
		xl = rtnewt_findtheta
		fl = f
	    end if
	    rtnewt_findtheta = 0.5d0*(xl+xh)
	  end do
	end if
!	do j=1,MAXIT
!		call funcd(rtnewt_findtheta,f,df,xrho,xtemp)
!		dx=f/df
!		rtnewt_findtheta=rtnewt_findtheta-dx
!		if ((x1-rtnewt_findtheta)*(rtnewt_findtheta-x2) < 0.0)&
!                     write(*,*) "rtnewt_findtheta: values jumped out of brackets"
!			!call nrerror('rtnewt_findtheta: values jumped out of brackets')
!		if (abs(dx) < xacc) RETURN
!	end do
	!call nrerror('rtnewt_findtheta exceeded maximum iterations')
        write(*,*) "rtnewt_findtheta exceeded maximum iterations"
	END FUNCTION rtnewt_findtheta
