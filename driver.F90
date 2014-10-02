! To Do (ordered by priority/dependency, with bottom depending on top):
!   * get the WARM tables' rootsolves to converge in especially difficult cases
!   * replace unit scalings with constants from eosmodule, check change in new tables
!   * replace table extents and bounds with those from eosmodule
!   * remove unneccesary arguments in functions and subroutines
!   * convert functions to subroutines
!   * move all helper subroutines to a separate module (driverhelpermodule)
!   * put driverhelpermodule in a separate *.F90 file
!   * include global constants (table type flags) in driverhelpermodule
!   * rewrite in cpp (just kidding)

program driver

  use eosmodule
  implicit none

  real*8 googol

  real*8 rtnewt_findYe
  real*8 rtnewt_findTYe
  real*8 rtnewt_findtheta
  
  real*8 xrho,xye,xtemp,tableymin
  real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
  real*8 xmu_mismatch,dxmu_mismatch
  real*8 xp_mismatch
  real*8 xdpderho,xdpdrhoe
  integer keytemp,keyerr

  ! for full eos call:
  real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  real*8 xxa,xxh,xxn,xxp

  integer i,j,k
  integer nr, nt, ny
  real*8 ltmax, ltmin, lymin, lymax, lrmin, lrmax
  real*8 lrmin_for_user_chooses_low_rho_bound, ltmin_for_user_chooses_low_temp_bound
  real*8 pressure_ratio_for_warm_tables
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
  real*8 rho_at_hmin,temp_at_hmin,ye_at_hmin,beta_ye_at_this_t
  integer FULL,BETA,COLD,WARM,BETAYE,COLDYE,WARMYE,MICRO,MUX,DERIVS,COLDMU,WARMT,WARMP,eostype
  integer HSHEN_2011,LS_220,GSHEN_NL3,GSHEN_FSU21,HEMPEL_SFHO,HEMPEL_SFHX,HEMPEL_DD2,eos
  logical USER_CHOOSES_BOUNDS
  logical USER_CHOOSES_LOW_RHO_BOUND
  logical USER_CHOOSES_LOW_TEMP_BOUND
  logical MONOTONIZE_MU_P
  logical MONOTONIZE_MU_N

  ! for safety checks on infinities
  googol=1.d100

  RTACC = 1.e-5      ! Relative accuracy for root of f(x). We use this in x and f both.
  keytemp = 1
  keyerr  = 2
  RhoScale = 5.e-5   ! this rescaling currently isn't used

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
  WARM = 4     ! internal energy in 1D table, with t!=tmin, ye=beta_equil
       	       !   This is different than a COLD table because t varies s.t. the thermal pressure is
	       !   a constant percentage of the total pressure at all densities. The thermal pressure
	       !   is defined as p(t)-p(t=tmin); so this table type will be different for different
	       !   choices of tmin.
               !   NOTE: WARM* solves take a while because they must perform a double-nested rootsolve.
  ! Tables of ye
  BETAYE = 5   ! Ye in 2D table, with ye=beta_equil
  COLDYE = 6   ! Ye in 1D table, with t=tmin, ye=beta_equil
  WARMYE = 7   ! Ye in 1D table, with t!=tmin, ye=beta_equil
  	       !   This follows the same t prescription as the WARM table.
  ! Tables of other thermodynamic potentials
  MICRO = 8    ! microphysical potentials in 3D table
  MUX = 9      ! chemical potentials in 3D table
  ! Testing Tables
  DERIVS = 10  ! thermodynamic derivs in 3D table
  COLDMU = 11  ! mu_mismatch in 1D table (varying Ye) to examine monotonicity of e,p,n potentials
  WARMT = 12   ! temperature in 1D table (varying rho) to examine the t rootsolve for WARM tables
  	       !   This follows the same t prescription as the WARM table.
  WARMP = 13   ! p_mismatch in 1D table (varying T) to examine the shape of the mismatch function
               !   for WARM tables.
  	       !   This follows the same t prescription as the WARM table.

  ! Define versions of tables:
  HSHEN_2011 = 1   ! HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5
  LS_220 = 2       ! LS220_234r_136t_50y_analmu_20091212_SVNr26.h5
  GSHEN_NL3 = 3    ! GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5
  GSHEN_FSU21 = 4  ! GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5
  HEMPEL_SFHO = 5  ! Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5
  HEMPEL_SFHX = 6  ! Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5
  HEMPEL_DD2 = 7   ! Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5

  ! ***** User-Chosen Parameters ************************************************************
  eostype = FULL
  eos = LS_220

  ! Choose the ratio of thermal pressure to total pressure to hold constant for WARM style tables.
  ! NOTE: if (eostype.ne.WARM).and.(eostype.ne.WARMYE).and.(eostype.ne.WARMT) then this parameter
  ! plays no role in the table generation.
  pressure_ratio_for_warm_tables = 1d-2

  ! Monotonize the proton and/or neutron chemical potential? Necessary for nonsmooth GShen tables.
  ! False if we want to use the chemical potential for p directly from the table.
  ! True is more computationally expensive, and also manipulates the tabulated mu_p and/or mu_n.
  ! I find that mu_p is the major problematic variable, and additionally monotonizing mu_n doesn't
  ! add much, except occasionally in the high-Ye region.
  MONOTONIZE_MU_P = .false.
  MONOTONIZE_MU_N = .false.

  USER_CHOOSES_LOW_RHO_BOUND = .true. ! true if user is to override table's low bound in r
  ! Choose low density bound different than table's intrinsic bound in r.
  !   This option is often used rather than USER_CHOOSES_BOUNDS because we want to honor
  !   the table's intrinsic bounds in all vars, except the minimum rho, where we do not need
  !   these lowest densities because we apply atmosphere treatment near rho=1e-10.
  !   if USER_CHOOSES_LOW_RHO_BOUND then user must supply low density bound here
  !   if .not.USER_CHOOSES_LOW_RHO_BOUND then this is never used
  lrmin_for_user_chooses_low_rho_bound = 5d0

  USER_CHOOSES_LOW_TEMP_BOUND = .false. ! true if user is to override table's low bound in t
  ! Choose low temperature bound different than table's intrinsic bound in t.
  !   This option is sometimes used rather than USER_CHOOSES_BOUNDS because we want to honor
  !   the table's intrinsic bounds in all vars, except the minimum temp, where we may not trust
  !   the lowest temperatures. This option is usually employed to make ColdTables by slicing
  !   the 3D table at a higher temperature than minimum.
  !   if USER_CHOOSES_LOW_TEMP_BOUND then user must supply low temperature bound here
  !   if .not.USER_CHOOSES_LOW_RHO_BOUND then this is never used
  ltmin_for_user_chooses_low_temp_bound = -0.7d0

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

  if (USER_CHOOSES_BOUNDS.and.USER_CHOOSES_LOW_TEMP_BOUND) then
    write(*,"(A,I1,A)") 'Error, user has specified a lower temp bound twice.'
  end if

  ! choose output resolution (nt and/or ny overwritten for some eostypes)
  nr = 345 ! use higher resolution for cold tables
  !nr = 250
  nt = 120
  ny = 100
  !ny = 2000 ! use higher resolution for COLDMU table

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

  ! reset ltmin to user-chosen bound if this option is used
  if (USER_CHOOSES_LOW_TEMP_BOUND) then
     ltmin = ltmin_for_user_chooses_low_temp_bound
  end if

  rmin = 10**lrmin
  rmax = 10**lrmax
  tmin = 10**ltmin
  tmax = 10**ltmax
  lymin = log10(ymin)
  lymax = log10(ymax)

  ! reset output resolutions for lower-dimensional tables
  if((eostype.eq.COLD).or.(eostype.eq.WARM).or.(eostype.eq.COLDYE).or.(eostype.eq.WARMYE).or.(eostype.eq.WARMT)) then
     nt = 1
     ny = 1
  end if
  if(eostype.eq.COLDMU) then
     nt = 1
     nr = 1
  end if
  if(eostype.eq.WARMP) then
     nr = 1
     ny = 1
  end if
  if((eostype.eq.BETA).or.(eostype.eq.BETAYE)) then
     ny = 1
  end if

  ! ***** Print Header: nice things for humans to know ***********************************************************
  if(eostype.eq.FULL) then
    write(*,"(A)") "# 3D full table: energy, pressure, sound speed"
  else if(eostype.eq.BETA) then
    write(*,"(A)") "# 2D beta-equilibrium table: energy, pressure, sound speed"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e = 0"
  else if(eostype.eq.COLD) then
    write(*,"(A)") "# 1D cold beta-equilibrium table: energy, pressure, sound speed"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e = 0"
    write(*,"(A,E15.6)") "# T sliced so that T = ", tmin
  else if(eostype.eq.WARM) then
    write(*,"(A)") "# 1D warm beta-equilibrium table: energy, pressure, sound speed"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e=0"
    write(*,"(A,E15.6,A)") "# T sliced so that P_thermal = ", pressure_ratio_for_warm_tables, "*P_total"
  else if(eostype.eq.BETAYE) then
    write(*,"(A)") "# 2D beta-equilibrium table: ye"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e = 0"
  else if(eostype.eq.COLDYE) then
    write(*,"(A)") "# 1D cold beta-equilibrium table: rho, ye"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e = 0"
    write(*,"(A,E15.6)") "# T sliced so that T = ", tmin
  else if(eostype.eq.WARMYE) then
    write(*,"(A)") "# 1D warm beta-equilibrium table: rho, ye"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e=0"
    write(*,"(A,E15.6,A)") "# T sliced so that P_thermal = ", pressure_ratio_for_warm_tables, "*P_total"
  else if(eostype.eq.MICRO) then
    write(*,"(A)") "# 3D full table: microphysical potentials"
  else if(eostype.eq.MUX) then
    write(*,"(A)") "# 3D full table: chemical potentials"
  else if(eostype.eq.DERIVS) then
    write(*,"(A)") "# 3D full table: thermodynamic derivatives"
    write(*,"(A)") "# FOR TESTING/DEBUGGING"
  else if(eostype.eq.COLDMU) then
    write(*,"(A)") "# 1D cold beta-equilibrium table: Ye, mu, dmu/dYe"
    write(*,"(A)") "# FOR TESTING/DEBUGGING"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e = 0"
    write(*,"(A,E15.6)") "# T sliced so that T = ", tmin
    write(*,"(A,E15.6)") "# rho sliced so that rho = ",rmin/rhounit
    write(*,"(A)") "# [1] Ye"
    write(*,"(A)") "# [2] mu=mu_n-mu_p-mu_e"
    write(*,"(A)") "# [3] dmu/dYe"
  else if(eostype.eq.WARMT) then
    write(*,"(A)") "# 1D warm beta-equilibrium table: rho, temp"
    write(*,"(A)") "# FOR TESTING/DEBUGGING"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e=0"
    write(*,"(A,E15.6,A)") "# T sliced so that P_thermal = ", pressure_ratio_for_warm_tables, "*P_total"
    write(*,"(A)") "# [1] rho"
    write(*,"(A)") "# [2] temp"
  else if(eostype.eq.WARMP) then
    write(*,"(A)") "# 1D warm beta-equilibrium table: temp, p_mismatch"
    write(*,"(A)") "# FOR TESTING/DEBUGGING"
    write(*,"(A)") "# Ye sliced so that mu_p-mu_n+mu_e=0"
    write(*,"(A,E15.6,A)") "# T sliced so that P_thermal = ", pressure_ratio_for_warm_tables, "*P_total"
    write(*,"(A)") "# [1] temp"
    write(*,"(A)") "# [2] p_mismatch"
  else
    write(*,"(A,I1,A)") 'Error, eostype ', eostype, ' is not recognized.'    
  end if

  ! ***** Print Header: warning about table manipulation ******************************************************
  if ((eostype.eq.BETA).or.(eostype.eq.COLD).or.(eostype.eq.WARM).or.(eostype.eq.BETAYE) &
      .or.(eostype.eq.COLDYE).or.(eostype.eq.WARMYE).or.(eostype.eq.COLDMU).or.(eostype.eq.WARMT) &
      .or.(eostype.eq.WARMP)) then
    if (MONOTONIZE_MU_P) then
      write(*,"(A)") "# Note: beta equilibrium Ye-curve smoothed by artificially enforcing monotonic mu_p(Ye) relation"
    end if
    if (MONOTONIZE_MU_N) then
      write(*,"(A)") "# Note: beta equilibrium Ye-curve smoothed by artificially enforcing monotonic mu_n(Ye) relation"
    end if
  end if

  ! ***** Print Header: SpEC classname ***********************************************************************
  if(eostype.eq.FULL) then
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.BETA) then
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.COLD) then
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.WARM) then
    write(*,"(A)") "EoSType = Tabulated"
  else if(eostype.eq.BETAYE) then
    write(*,"(A)") "EoSType = HotBetaYe"
  else if(eostype.eq.COLDYE) then
    write(*,"(A)") "EoSType = ColdBetaYe"
  else if(eostype.eq.WARMYE) then
    write(*,"(A)") "EoSType = ColdBetaYe"
  else if(eostype.eq.MICRO) then
    write(*,"(A)") "EoSType = Microphysics"
  else if(eostype.eq.MUX) then
    write(*,"(A)") "EosType = ChemicalPotentials"
  else if((eostype.eq.DERIVS).or.(eostype.eq.COLDMU).or.(eostype.eq.WARMT).or.(eostype.eq.WARMP)) then
    write(*,"(A)") "# EosType = this isn't a SpEC-able table, gnuplot it"
  else
    write(*,"(A,I1,A)") 'Error, eostype ', eostype, ' is not recognized.'    
  end if

  ! ***** Print Header: extents ****************************************************************************
  if(eostype.eq.BETAYE) then
    write(*,"(A,I6,A,A,I6)") "Nrho =",nr,"     ","NT =",nt
  else if((eostype.eq.COLDYE).or.(eostype.eq.WARMYE)) then
    write(*,"(A,I6)") "Nrho =",nr
  else if(eostype.eq.DERIVS) then
    write(*,"(A,I6,A,A,I6,A,A,I6)") "# Nrho =",nr,"     ","Nye =",ny,"     ","NT =",nt
  else if(eostype.eq.COLDMU) then
    write(*,"(A,I6)") "# Nye =",ny
  else if(eostype.eq.WARMT) then
    write(*,"(A,I6)") "# Nrho =",nr
  else if(eostype.eq.WARMP) then
    write(*,"(A,I6)") "# NT =",nt
  else  
    write(*,"(A,I6,A,A,I6,A,A,I6)") "Nrho =",nr,"     ","Nye =",ny,"     ","NT =",nt
  end if

  ! ***** Print Header: bounds ****************************************************************************
  ! Table with rho bounds
  if((eostype.ne.COLDMU).and.(eostype.ne.WARMP)) then
    write(*,"(A, E15.6, A, A, E15.6)") "RhoMin =",rmin/rhounit,"      ",&
    		 	      	       "RhoMax =",rmax/rhounit
  end if
  ! Tables with ye bounds
  if((eostype.eq.FULL).or.(eostype.eq.MICRO).or.(eostype.eq.MUX) &
     .or.(eostype.eq.DERIVS).or.(eostype.eq.COLDMU)) then
    write(*,"(A, E15.6, A, A, E15.6)") "YeMin =",ymin,"      ",&
    		 	      	       "YeMax =",ymax
  end if
  ! Tables with T bounds
  if((eostype.eq.FULL).or.(eostype.eq.BETA).or.(eostype.eq.BETAYE) &
     .or.(eostype.eq.MICRO).or.(eostype.eq.MUX).or.(eostype.eq.DERIVS).or.(eostype.eq.WARMP)) then
    write(*,"(A, E15.6, A, A, E15.6)") "Tmin =",tmin,"      ",&
    		 	      	       "Tmax =",tmax
  end if

  ! ***** Print Header: other metadata *********************************************************************  
  kappa = 100.d0        ! for our output of gamma s.t. pressure=kappa*rho^gamma
  HeatCapacity = 1.6d-3 ! currently used as an extrapolation parameter for temperature
  GammaTh = 5.d0/3.d0   ! same ^
  gtfac = (GammaTh-1.d0)/GammaTh
  if((eostype.eq.FULL).or.(eostype.eq.BETA).or.(eostype.eq.COLD) &
       .or.(eostype.eq.WARM)) then
     write(*,"(A, E15.6)") "HeatCapacity =",HeatCapacity,&
          "GammaTh =",GammaTh,&
          "Kappa =",kappa
  end if

  ! ***** Print Header: spacing ***************************************************************************  
  ! Tables with rho bounds
  if((eostype.ne.COLDMU).and.(eostype.ne.WARMP)) then
    write(*,"(A)") "RhoSpacing = Log"
  end if
  ! Tables with ye bounds
  if((eostype.eq.FULL).or.(eostype.eq.MICRO).or.(eostype.eq.MUX) &
     .or.(eostype.eq.DERIVS).or.(eostype.eq.COLDMU)) then
    write(*,"(A)") "YeSpacing = Linear"
  end if
  ! Tables with T bounds
  if((eostype.eq.FULL).or.(eostype.eq.BETA).or.(eostype.eq.BETAYE) &
     .or.(eostype.eq.MICRO).or.(eostype.eq.MUX).or.(eostype.eq.DERIVS).or.(eostype.eq.WARMP)) then
    write(*,"(A)") "TSpacing = Log"
  end if

  ! ***** Print Table ************************************************************************
  ! Here we print the input table used in a SpEC Hydro evolution.
  !   To be read in using the EquationOfState::Tabulated class.
  !   For the standard tables -- FULL (3D), BETA (2D), COLD(1D), WARM(1D) -- as of May 2013,
  !   the output colums are:
  ! [1] epsilon, specific internal energy
  ! [2] gamma,   so that P=kappa*rho^gamma, where kappa is the constant defined above
  ! [3] cs,      relativistic adiabatic sound speed

  ! In the current implementation hmin vars are not used; they were used once to find hmin in the table.  
  hmin = 1.d0
  rho_at_hmin = rmin
  temp_at_hmin = tmin
  ye_at_hmin = ymin
  ! T loop
  do i=0, nt-1
    if(nt.eq.1) then
      xtemp = tmin
    else
      !xtemp = tmin + i*(tmax-tmin)/(nt-1)         ! linear
      xtemp = 10**(ltmin + i*(ltmax-ltmin)/(nt-1)) ! logarithmic
    end if
    ! Ye loop
    do j=0, ny-1
      if(ny.eq.1) then
        xye = ymin
      else
        xye = ymin + j*(ymax-ymin)/(ny-1)           ! linear
        !xye = 10**(lymin + j*(lymax-lymin)/(ny-1)) ! logarithmic
      end if
      ! rho loop
      do k=0, nr-1
        if(nr.eq.1) then
	  xrho = rmin
	else
          xrho = 10**(lrmin + k*(lrmax-lrmin)/(nr-1)) ! logarithmic
	end if

	! set xye to satisfy beta_equil
	! We don't have to do this for WARM* tables, because that's done within rtnewt_findTYe
        if((eostype.eq.BETA).or.(eostype.eq.COLD) &
	   .or.(eostype.eq.BETAYE).or.(eostype.eq.COLDYE) &
	   .or.(eostype.eq.COLDMU)) then
          xye = rtnewt_findYe(ymin,ymax,RTACC,xrho,xtemp,tableymin,MONOTONIZE_MU_P,MONOTONIZE_MU_N)
        end if

	! now query the table at this rho,T,Ye
        call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)

        ! sometimes xcs2 is Infinity (at very low densities and low temps),
        !   where it should be vanishing)
        if(xcs2.gt.googol) then
           xcs2=0.d0
        end if

	! find the microphysical potential theta
        if(eostype.eq.MICRO) then
	  xrhothet = 3.e12
          theta = rtnewt_findtheta(tmin,tmax,RTACC,xrhothet,xye,xent)
        end if

        ! find the mismatch in e,p,n chemical potentials for COLDMU table
	if(eostype.eq.COLDMU) then
          call mu_mismatch(xye,xmu_mismatch,dxmu_mismatch,xrho,xtemp,tableymin,MONOTONIZE_MU_P,MONOTONIZE_MU_N)
	end if

        ! find the mismatch in relative thermal pe and target relative thermal p for WARMP table
	if(eostype.eq.WARMP) then
          call p_mismatch(xtemp,xp_mismatch,xrho,xye,tableymin,ymax,RTACC,tmin,tableymin, &
               pressure_ratio_for_warm_tables,MONOTONIZE_MU_P,MONOTONIZE_MU_N)
	end if

	! Find the t (at ye=beta_equilibrium) which keeps the thermal pressure at a constant
	!   fraction of the total pressure.
	if((eostype.eq.WARM).or.(eostype.eq.WARMYE).or.(eostype.eq.WARMT)) then
	  beta_ye_at_this_t = ymin
          ! silly hack (putting rtnewt_findTYe's return value into xtemp_p instead of xtemp) necessary
          !   for now, until I convert this function to a subroutine. We don't use xtemp_p for these tables.
	  xtemp_p = rtnewt_findTYe(xtemp,beta_ye_at_this_t,tmin,tmax,RTACC,tableymin,ymax,RTACC, &
				 xrho,tmin,tableymin,pressure_ratio_for_warm_tables,MONOTONIZE_MU_P,MONOTONIZE_MU_N)
	  xye = beta_ye_at_this_t
          call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
               xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
               xmuhat,keytemp,keyerr)
	end if

	! take centered derivs
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

	! convert to code units G=c=Msun=1
        xrho = xrho/rhounit
        xenr = xenr/eomunit
        xcs2 = xcs2/c2unit
        xprs = xprs/punit
	xent = xent*9.56565348d17/eomunit

        GammaP = log(xprs/kappa)/log(xrho)
        EoverRho = xenr

	! scale eps to zero
        if(xxn<1.e-20) xxn = 0.d0
        if(xxp<1.e-20) xxp = 0.d0
        if(xxa<1.e-20) xxa = 0.d0
        if(xxh<1.e-20) xxh = 0.d0

	! keep track of minimum enthalpy
	if(1.0+xenr+xprs/xrho.lt.hmin) then
	  hmin = 1.0+xenr+xprs/xrho
	  rho_at_hmin = xrho
	  temp_at_hmin = xtemp
	  ye_at_hmin = xye
	end if

	! Print the potentials/primitives of interest
        if((eostype.eq.FULL).or.(eostype.eq.BETA) &
	   .or.(eostype.eq.COLD).or.(eostype.eq.WARM)) then
	  cs = 0.d0
	  if(xcs2.gt.0.d0) cs = sqrt(xcs2/(1.0+xenr+xprs/xrho))
          write(*,"(1P10E15.8)") EoverRho,GammaP,cs
	else if(eostype.eq.BETAYE) then
	  write(*,"(1P10E15.6)") xye
	else if(eostype.eq.COLDYE) then
	  write(*,"(1P10E15.6)") xrho, xye
	else if(eostype.eq.WARMYE) then
	  write(*,"(1P10E15.6)") xrho, xye
        else if(eostype.eq.MICRO) then
          write(*,"(1P10E15.6)") xent,xxn,xxp,xxa,xxh,xabar,xzbar,theta	
	else if(eostype.eq.MUX) then
	  write(*,"(1P10E15.6)") xmu_e,xmu_p-xmu_n,xxn,xxp,xxh,xabar,xzbar
	else if(eostype.eq.DERIVS) then
	  write(*,"(1P10E15.6)") pgam,dPds,dPdY,xent
	else if(eostype.eq.COLDMU) then
	  write(*,"(1P10E15.6)") xye, xmu_mismatch, dxmu_mismatch
	else if(eostype.eq.WARMT) then
	  write(*,"(1P10E15.6)") xrho, xtemp
	else if(eostype.eq.WARMP) then
	  write(*,"(1P10E15.6)") xtemp, xp_mismatch
	else
          write(*,"(A,I1,A)") 'Error, eostype ', eostype, ' is not recognized.'
        end if
      end do ! END for r
    end do ! END for y
  end do ! END for t
end program driver

! ****************************************************************************************

! Calculate a monotonic version of the mu_p curve.
! In some EOSs (notably GShenNL3 and GShenFSU21) the mu_p curve at a given rho,Temp
! is nonmonotonic in Ye, giving us poor root-solves with the newton-raphson method.
! We resolve this by scanning the mu_p(Ye) curve at each point in rho,Temp, from high Ye
! to low; we keep mu_p constant if the tabulated value is nonmonotonic, otherwise, we use
! the tabulated value.
! Since we only have to do this to a few tables, this is coded up *incredibly*
! inefficiently: for each root solve call at the same rho,Temp we recalculate the entire
! mu_p(Ye) curve. We only bother to do this at low-ish electron fractions,
! where the root-solver spends most of its iterations.
! return value:              xmu_p
! table abscissa:            xye,xrho,xtemp
! table low bound in Ye:     tableymin
! highest Ye to monotonize:  max_monotonized_ye
SUBROUTINE monotonized_xmu_p(xmu_p,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xrho,xtemp,xye,tableymin,max_monotonized_ye
  REAL*8, INTENT(OUT) :: xmu_p
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmuhat
  integer keytemp,keyerr
  real*8 :: xye_j,xye_jj                  ! j is the lower point at (j), jj is the upper point at (j+1)
  real*8 :: xmu_p_j,xmu_p_jj,xmu_p_temp   ! xmu_p_temp is our 'peak ahead' variable
  real*8 :: slope                         ! slope of mu_p,ye curve for linear interpolation
  integer ny, j
  keytemp = 1
  keyerr = 0

  if(xye<max_monotonized_ye) then
    ! Calculate the monotonized curve
    ! Since most input EOS tables have about 50 points in Ye between 0 and 0.5, we'll make
    ! a grid of 50 points between tableymin and max_monotonized_ye.
    ! From this grid, we'll find the points that bound xye, and linearly interpolate.
    ny=50
    ! set mu_p and ye for upper point
    xye_jj = max_monotonized_ye
    call nuc_eos_full(xrho,xtemp,xye_jj,xenr,xprs,xent,xcs2,xdedt,&
      xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p_jj,&
      xmuhat,keytemp,keyerr)
    xmu_p_j = xmu_p_jj ! initialize the mu_p(ye) slope to zero
    ! start the loop at the penultimate point
    do j=1, ny
      xye_j = max_monotonized_ye - j*(max_monotonized_ye-tableymin)/(ny-1)
      call nuc_eos_full(xrho,xtemp,xye_j+1d-12,xenr,xprs,xent,xcs2,xdedt,&
        xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p_temp,&
        xmuhat,keytemp,keyerr) ! add a little eps to xye_j so we don't go below tableymin

      ! enforce monotonic decrease in xmu_p_j
      if (xmu_p_temp<xmu_p_j) then
        xmu_p_j = xmu_p_temp
      end if

      if (xye_j<xye) exit ! exit now, because we have the points that are bracketing xye

      ! if xye_j>xye, then shift the value in xmu_p_j to xmu_p_jj, (and same for ye) and continue loop
      xmu_p_jj = xmu_p_j
      xye_jj = xye_j
    end do

    ! interpolate to xye between xmu_p_j and xmu_p_jj
    slope = (xmu_p_jj-xmu_p_j) / (xye_jj-xye_j)
    xmu_p = slope*(xye-xye_j) + xmu_p_j

  else
    ! if xye>max_monotonized_ye, just use the raw table value
    call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
      xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
      xmuhat,keytemp,keyerr)
  end if

END SUBROUTINE monotonized_xmu_p

! ****************************************************************************************

! Calculate a monotonic version of the mu_n curve.
! In some EOSs (notably GShenNL3 and GShenFSU21) the mu_n curve at a given rho,Temp
! is nonmonotonic in Ye, giving us poor root-solves with the newton-raphson method.
! We resolve this by scanning the mu_n(Ye) curve at each point in rho,Temp, from high Ye
! to low; we keep mu_n constant if the tabulated value is nonmonotonic, otherwise, we use
! the tabulated value.
! Since we only have to do this to a few tables, this is coded up *incredibly*
! inefficiently: for each root solve call at the same rho,Temp we recalculate the entire
! mu_n(Ye) curve. We only bother to do this at low-ish electron fractions,
! where the root-solver spends most of its iterations.
! return value:              xmu_n
! table abscissa:            xye,xrho,xtemp
! table low bound in Ye:     tableymin
! highest Ye to monotonize:  max_monotonized_ye
SUBROUTINE monotonized_xmu_n(xmu_n,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xrho,xtemp,xye,tableymin,max_monotonized_ye
  REAL*8, INTENT(OUT) :: xmu_n
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: xye_j,xye_jj                  ! j is the lower point at (j), jj is the upper point at (j+1)
  real*8 :: xmu_n_j,xmu_n_jj,xmu_n_temp   ! xmu_n_temp is our 'peak ahead' variable
  real*8 :: slope                         ! slope of mu_n,ye curve for linear interpolation
  integer ny, j
  keytemp = 1
  keyerr = 0

  if(xye<max_monotonized_ye) then
    ! Calculate the monotonized curve
    ! Since most input EOS tables have about 50 points in Ye between 0 and 0.5, we'll make
    ! a grid of 50 points between tableymin and max_monotonized_ye.
    ! From this grid, we'll find the points that bound xye, and linearly interpolate.
    ny=50
    ! set mu_n and ye for upper point
    xye_jj = max_monotonized_ye
    call nuc_eos_full(xrho,xtemp,xye_jj,xenr,xprs,xent,xcs2,xdedt,&
      xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n_jj,xmu_p,&
      xmuhat,keytemp,keyerr)
    xmu_n_j = xmu_n_jj ! initialize the mu_n(ye) slope to zero
    ! start the loop at the penultimate point
    do j=1, ny
      xye_j = max_monotonized_ye - j*(max_monotonized_ye-tableymin)/(ny-1)
      call nuc_eos_full(xrho,xtemp,xye_j+1d-12,xenr,xprs,xent,xcs2,xdedt,&
        xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n_temp,xmu_p,&
        xmuhat,keytemp,keyerr) ! add a little eps to xye_j so we don't go below tableymin

      ! enforce monotonic increase in xmu_n_j
      if (xmu_n_temp>xmu_n_j) then
        xmu_n_j = xmu_n_temp
      end if

      if (xye_j<xye) exit ! exit now, because we have the points that are bracketing xye

      ! if xye_j>xye, then shift the value in xmu_n_j to xmu_n_jj, (and same for ye) and continue loop
      xmu_n_jj = xmu_n_j
      xye_jj = xye_j
    end do

    ! interpolate to xye between xmu_n_j and xmu_n_jj
    slope = (xmu_n_jj-xmu_n_j) / (xye_jj-xye_j)
    xmu_n = slope*(xye-xye_j) + xmu_n_j

  else
    ! if xye>max_monotonized_ye, just use the raw table value
    call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
      xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
      xmuhat,keytemp,keyerr)
  end if

END SUBROUTINE monotonized_xmu_n

! ****************************************************************************************

! Calculate mismatch of chemical potentials: fval=mu_n-mu_e-mu_p.
! return values:           fval, fderiv
! table abscissa:          xye, xrho, xtemp
! table low bound in Ye:   tableymin
! switch to care for mu_p: monotonize_mu_p,monotonize_mu_n
SUBROUTINE mu_mismatch(xye,fval,fderiv,xrho,xtemp,tableymin,monotonize_mu_p,monotonize_mu_n)
  !USE nrtype
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xye,xrho,xtemp,tableymin
  LOGICAL, INTENT(IN) :: monotonize_mu_p,monotonize_mu_n
  REAL*8, INTENT(OUT) :: fval,fderiv
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: ye_p, ye_m, EPS, f_p, f_m, deriv_crit, ymin, lymin, max_monotonized_ye
  EPS = 1.d-7
  keytemp = 1
  keyerr  = 0

  ! lookup chemical potentials
  call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)

  ! get the monotonized mu_p if that option flagged
  max_monotonized_ye=0.52d0
  if(monotonize_mu_p) then
    call monotonized_xmu_p(xmu_p,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if
  if(monotonize_mu_n) then
    call monotonized_xmu_n(xmu_n,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if

  fval = xmu_n - xmu_p - xmu_e

  ! set lower values for symmetric deriv
  ye_m = xye - EPS
  if(ye_m.lt.tableymin) ye_m = tableymin
  call nuc_eos_full(xrho,xtemp,ye_m,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  if(monotonize_mu_p) then
    call monotonized_xmu_p(xmu_p,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if
  if(monotonize_mu_n) then
    call monotonized_xmu_n(xmu_n,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if
  f_m = xmu_n - xmu_p - xmu_e

  ! set upper values for symmetric deriv
  ye_p = ye_m + 2.d0*EPS
  call nuc_eos_full(xrho,xtemp,ye_p,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  if(monotonize_mu_p) then
    call monotonized_xmu_p(xmu_p,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if
  if(monotonize_mu_n) then
    call monotonized_xmu_n(xmu_n,xrho,xtemp,xye,tableymin,max_monotonized_ye)
  end if
  f_p = xmu_n - xmu_p - xmu_e

  fderiv = (f_p - f_m)/(2.d0*EPS)

  ! keep next guess within ye-bounds, by adjusting the derivative at the boundaries
  deriv_crit = 1.01d0*fval/(xye-tableymin)
  if(abs(fderiv).lt.abs(deriv_crit)) then
    fderiv = deriv_crit
  end if
END SUBROUTINE mu_mismatch

! ****************************************************************************************

! Bisection root finding (name rtnewt is legacy) for mu_mismatch function
! function to search for root:    mu_mismatch
! table abscissa:                 xrho, xtemp
! table bounds in Ye:             x1, x2
! table low bound in Ye:          tableymin (only used to pass to mu_mismatch)
! requested mu_mismatch accuracy: xacc
! switch to care for mu_p,mu_n:   monotonize_mu_p, monotonize_mu_n
FUNCTION rtnewt_findYe(x1,x2,xacc,xrho,xtemp,tableymin,monotonize_mu_p,monotonize_mu_n)
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: x1,x2,xacc,xrho,xtemp,tableymin
  LOGICAL, INTENT(IN) :: monotonize_mu_p,monotonize_mu_n
  REAL*8 :: rtnewt_findYe
  INTEGER, PARAMETER :: MAXIT=40 ! 40 iterations of bisection gives an accuracy in x of 1e-12
  INTEGER :: j
  REAL*8 :: df,dx,f,fl,fh,dfl,dfh, xl,xh
  
  if(x1<x2) then
     xl = x1
     xh = x2
  else
     xl = x2
     xh = x1
  end if
  call mu_mismatch(xl,fl,dfl,xrho,xtemp,tableymin,monotonize_mu_p,monotonize_mu_n)
  call mu_mismatch(xh,fh,dfh,xrho,xtemp,tableymin,monotonize_mu_p,monotonize_mu_n)
  
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
        call mu_mismatch(rtnewt_findYe,f,df,xrho,xtemp,tableymin,monotonize_mu_p,monotonize_mu_n)
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
  write(*,*) "rtnewt_findYe exceeded maximum iterations, which is expected when the"
  write(*,*) "mu_mismatch function varies steeply with Ye. It is okay, the Ye value"
  write(*,*) "has been found with bisection to sufficient accuracy."
END FUNCTION rtnewt_findYe

! ****************************************************************************************

! Calculate mismatch of thermal pressure fraction and target: fval=(1-R)P(T)-P(Tmin).
!   The value of this function is scaled by a large factor (~1e30) toward unity,
!   because pressure in cgs is large.
! return values:           fval,this_ye
! table abscissa:          xrho, this_ye (beta-equilibrium ye, constrained variable), t
! accuracy of Ye rootfind: yacc
! table low bound in T,Ye: tabletmin, tableymin
! bounds of Ye root find:  ymin, ymax
! targe pressure ratio:    R (thermal pressure should be this fraction of total pressure)
! switch to care for mu_p: monotonize_mu_p, monotonize_mu_n
SUBROUTINE p_mismatch(t,fval,xrho,this_ye,ymin,ymax,yacc,tabletmin,tableymin,R,monotonize_mu_p,monotonize_mu_n)
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: t,xrho,tabletmin,tableymin,ymin,ymax,yacc,R
  LOGICAL, INTENT(IN) :: monotonize_mu_p,monotonize_mu_n
  REAL*8, INTENT(OUT) :: fval,this_ye
  REAL*8 :: rtnewt_findYe
  real*8 :: xye,xenr,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: p_at_t, p_at_tmin
  REAL*8 :: EPS    = 1.d-7        ! this is the step size for a numerical derivative
  REAL*8 :: PSCALE = 1.d30        ! this is the rescaling factor to move this function toward unity
  keytemp = 1
  keyerr  = 0

  ! get beta-equilibrium ye at this rho and t
  this_ye = rtnewt_findYe(ymin,ymax,yacc,xrho,t,tableymin,monotonize_mu_p,monotonize_mu_n)

  ! get the pressure at this rho,ye,t
  call nuc_eos_full(xrho,t,this_ye,xenr,p_at_t,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  ! get the pressure at tmin
  call nuc_eos_full(xrho,tabletmin,this_ye,xenr,p_at_tmin,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)

  fval = ((1-R)*p_at_t - p_at_tmin)/PSCALE

END SUBROUTINE p_mismatch

! ****************************************************************************************

! Bisection (name rtnewt is legacy) root finding for p_mismatch function, this is a 2D function,
!   defined over T and Ye. But we're iterating 2 1D rootfinds: Ye inside of T.
!   This function, different from rtnewt_findYe, is written to retrieve the return values (T,Ye)
!   from the argument list, not from the function value itself. Although, the function itself
!   does return the value of T that make the P mismatch vanish.
! return values, roots:           this_t, this_ye
! function to search for root:    p_mismatch
! table abscissa:                 xrho
! table bounds in T:              t1, t2
! table bounds in Ye:             y1, y2
! table low bound in T:           tabletmin (only used to pass to p_mismatch)
! table low bound in Ye:          tableymin (only used to pass to mu_mismatch)
! requested accuracy:             tacc, yacc
! target pressure ratio:          R (thermal pressure should be this fraction of total pressure)
! switch to care for mu_p,mu_n:   monotonize_mu_p, monotonize_mu_n
FUNCTION rtnewt_findTYe(this_t,this_ye,t1,t2,tacc,y1,y2,yacc,xrho,tabletmin,tableymin,R,monotonize_mu_p,monotonize_mu_n)
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: t1,t2,tacc,y1,y2,yacc,xrho,tabletmin,tableymin,R
  LOGICAL, INTENT(IN) :: monotonize_mu_p,monotonize_mu_n
  REAL*8, INTENT(OUT) :: this_t,this_ye
  REAL*8 :: rtnewt_findTYe
  INTEGER, PARAMETER :: MAXIT=40 ! 40 iterations of bisection gives an accuracy in x of 1e-12
  INTEGER :: j
  REAL*8 :: df,dt,dy,f,fl,fh,dfl,dfh,tl,th,yl,yh
  REAL*8 :: this_ye_h,this_ye_l
  
  if(t1<t2) then
     tl = t1
     th = t2
  else
     tl = t2
     th = t1
  end if
  
  ! fetch the initial bracket
  call p_mismatch(tl,fl,xrho,this_ye_l,y1,y2,yacc,tabletmin,tableymin,R,monotonize_mu_p,monotonize_mu_n)
  call p_mismatch(th,fh,xrho,this_ye_h,y1,y2,yacc,tabletmin,tableymin,R,monotonize_mu_p,monotonize_mu_n)
  
  ! Check if one of these limits is good enough already
  if(abs(fl)<tacc) then
     this_t = tl
     this_ye = this_ye_l
     RETURN
  end if
  if(abs(fh)<tacc) then
     this_t = th
     this_ye = this_ye_h
     RETURN
  end if
  
  if(fl*fh>0.d0) then
     this_t = 0              !
     rtnewt_findTYe = this_t ! I'm not sure if we have to set the value of the function or not
     write(*,"(A, E15.6, A, E15.6)") "rtnewt_findTYe failed to find an initial p_mismatch bracket in T between ", &
          tl, ", and ", th
     RETURN
  else 
     this_t = 0.5d0*(tl+th)               ! make initial bisection
     do j=1,MAXIT
        call p_mismatch(this_t,f,xrho,this_ye,y1,y2,yacc,tabletmin,tableymin,R,monotonize_mu_p,monotonize_mu_n)
        if(abs(f)<tacc) then
           rtnewt_findTYe = this_t
           RETURN
        end if
        if(f*fl<0.d0) then
           th = this_t
           fh = f
        else
           tl = this_t
           fl = f
        end if
        this_t = 0.5d0*(tl+th)             ! make next bisection
        rtnewt_findTYe = this_t            ! in case we exceed max it, this assignment necessary
     end do
  end if
  write(*,*) "rtnewt_findTYe exceeded maximum iterations, which is expected when the"
  write(*,*) "pressure mismatch function varies dramatically with temperature. It is okay,"
  write(*,*) "the temperature has been found to sufficient accuracy with bisection."
  
END FUNCTION rtnewt_findTYe

! ****************************************************************************************

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

  call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)

  fval = xent - xentGoal
  tmin = 0.1d0
  deriv_crit = 1.01d0*fval/(xtemp-tmin)
  if(abs(fderiv).lt.abs(deriv_crit)) then
    fderiv = deriv_crit
  end if

END SUBROUTINE ent_mismatch

! ****************************************************************************************

! Bisection root find (name rtnewt is legacy) over the Theta function
FUNCTION rtnewt_findtheta(x1,x2,xacc,xrho,xye,xent)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xye,xent
  REAL*8 :: rtnewt_findtheta
  INTEGER, PARAMETER :: MAXIT=40
  INTEGER :: j
  REAL*8 :: df,dx,f,fl,fh,dfl,dfh,xl,xh

  !	rtnewt_findtheta=0.5d0*(x1+x2)
  !	call ent_mismatch(rtnewt_findtheta,f,df,xrho,xye,xent)
  if(x1<x2) then
     xl = x1
     xh = x2
  else
     xl = x2
     xh = x1
  end if
  call ent_mismatch(xl,fl,dfl,xrho,xye,xent)
  call ent_mismatch(xh,fh,dfh,xrho,xye,xent)
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
        call ent_mismatch(rtnewt_findtheta,f,df,xrho,xye,xent)
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

  write(*,*) "rtnewt_findtheta exceeded maximum iterations"

END FUNCTION rtnewt_findtheta
