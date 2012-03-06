program driver

  use eosmodule
  implicit none

  interface
	FUNCTION rtnewt_findYe(funcd,x1,x2,xacc,xrho,xtemp)
	!USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: x1,x2,xacc,xrho,xtemp
	REAL*8 :: rtnewt_findYe
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv,xrho,xtemp)
		!USE nrtype
		IMPLICIT NONE
		REAL*8, INTENT(IN)  :: x,xrho,xtemp
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
	SUBROUTINE mu_mismatch(x,fval,fderiv,xrho,xtemp)
	!USE nrtype
	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: x,xrho,xtemp
	REAL*8, INTENT(OUT) :: fval,fderiv
	END SUBROUTINE mu_mismatch
	SUBROUTINE ent_mismatch(x,fval,fderiv,xye,xentGoal)
	!USE nrtype
	IMPLICIT NONE
	REAL*8, INTENT(IN)  :: x,xye,xentGoal
	REAL*8, INTENT(OUT) :: fval,fderiv
	END SUBROUTINE ent_mismatch
  end interface
  real*8 xrho,xye,xtemp,xs
  real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
  real*8 xdpderho,xdpdrhoe
  integer keytemp,keyerr

  ! for full eos call:
  real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  real*8 xxa,xxh,xxn,xxp

  integer i,j,k
  integer nr, nt, ny
  real*8 ltmax, ltmin, lymin, lymax, lrmin, lrmax
  real*8 tmax, tmin, ymin, ymax, rmin, rmax, tminout, tmaxout
  real*8 yeavg
  real*8 xunit,tunit,munit,rhounit,punit,eomunit,c2unit,entunit
  real*8 GammaP, EoverRho
  real*8 RhoScale, HeatCapacity, GammaTh, kappa, gtfac
  real*8 BARYMASS, RTACC
  real*8 r_tmin, r_tmax
  real*8 theta, xrhothet, xentthet
  real*8 pgam,delta_t,delta_y,xtemp_p,xtemp_pp,xye_p
  real*8 xprs_tp,xent_tp,xenr_tp,xprs_tpp,xent_tpp,xenr_tpp
  real*8 xprs_yp,xent_yp,dPdT,dsdT,dPdY,dsdY,dPdYs,dPds
  real*8 hmin,cs
  real*8 ent_min,ent_max, lent_min,lent_max

  integer EOSTABLE,TEMP_TO_ENT,ENT_TO_TEMP,MICRO,MUX,output
  integer LS,SHEN,eos
  EOSTABLE    = 1
  TEMP_TO_ENT = 2
  ENT_TO_TEMP = 3
  MICRO       = 4
  MUX         = 5
  LS = 1
  SHEN = 2

  output = ENT_TO_TEMP
  eos = LS 

  RTACC = 1.e-5
  keytemp = 1
  keyerr  = 0

  nr = 250
  nt = 120
  ny = 100

  ltmin = -1.0d0
  ltmax = 1.95d0
  lymax = -0.3d0
  lrmin = 8d0
  lrmax = 15.d0
  tmin = 10**ltmin
  tmax = 10**ltmax
  if(eos.eq.SHEN) then
    ymin = 0.02d0
  else
    ymin = 0.035d0
  end if
  ymax = 10**lymax
  rmin = 10**lrmin
  rmax = 10**lrmax

  xunit = 1.4766d5
  tunit = 4.9255d-6
  munit = 1.989d33
  BARYMASS = 1.659d-24
  rhounit = munit/xunit**3
  punit = munit/(xunit*tunit**2)
  c2unit = (xunit/tunit)**2
  eomunit = c2unit
  entunit = eomunit/9.56565348d17

  RhoScale = 5.e-5

  if(eos.eq.SHEN) then
    call readtable("myhshen2010_220r_180t_63y_analmu_20100831.h5")
  else
    call readtable("LS220_450r_270t_50y_062211.h5")
  end if

  if(output.eq.MICRO) then
    write(*,*) "EoSType = Microphysics"
  else if(output.eq.MUX) then
    write(*,*) "EoSType = ChemicalPotentials"
  else if(output.eq.EOSTABLE) then
    write(*,*) "EoSType = Tabulated"
  else if(output.eq.TEMP_TO_ENT) then
    write(*,*) "EoSType = Temperature"
  else if(output.eq.ENT_TO_TEMP) then
    write(*,*) "EoSType = Temperature"
  end if

! find ent_min, ent_max
  call nuc_eos_full(rmax,tmin,ymin,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)
  ent_min = xent*sqrt(rmax/rhounit)
  call nuc_eos_full(rmax,tmin,ymax,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)
  if(xent*sqrt(rmax/rhounit)<ent_min) then
    ent_min = xent*sqrt(rmax/rhounit)
  end if
  call nuc_eos_full(rmin,tmax,ymin,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)
  ent_max = xent*sqrt(rmin/rhounit)
  call nuc_eos_full(rmin,tmax,ymax,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)
  if(xent*sqrt(rmin/rhounit)>ent_max) then
    ent_max = xent*sqrt(rmin/rhounit)
  end if

  lent_min = log10(ent_min)
  lent_max = log10(ent_max)

  write(6,"(A6,I4,A5,I4,A6,I4)") "Nrho=",nr,"Nye=",ny,"NS=",nt

  write(6,"(A8, E15.6, A11, E15.6)") "RhoMin=",rmin/rhounit, &
        "RhoMax=",rmax/rhounit
  write(6,"(A7, E15.6, A10, E15.6)") "YeMin=",ymin, &
        "YeMax=",ymax
  if(output.eq.TEMP_TO_ENT) then
    write(6,"(A6, E15.6, A9, E15.6)") "Tmin=",tmin, &
          "Tmax=",tmax
  else
    write(6,"(A6, E15.6, A9, E15.6)") "Smin=",ent_min/entunit, &
          "Smax=",ent_max/entunit
  end if

  if(output.eq.EOSTABLE) then
    kappa = 100.d0
    HeatCapacity = 1.6d-3
    GammaTh = 5.d0/3.d0
    gtfac = (GammaTh-1.d0)/GammaTh
    write(6,"(A13, E15.6)") "HeatCapacity=",HeatCapacity," GammaTh=",GammaTh
    write(6,"(A10, E15.6)") "Kappa=",kappa
  end if

  write(*,*) "RhoSpacing = Log"
  write(*,*) "YeSpacing = Linear"
  write(*,*) "SSpacing = Log"

  hmin = 1.d0
  do i=0, nt-1
     if(output.eq.TEMP_TO_ENT) then
	xtemp = 10**(ltmin + i*(ltmax-ltmin)/(nt-1))
     else
        xs = 10**(lent_min + i*(lent_max-lent_min)/(nt-1))
     end if
   do j=0, ny-1
      xye = ymin + j*(ymax-ymin)/(ny-1)
      do k=0, nr-1
        xrho = 10**(lrmin + k*(lrmax-lrmin)/(nr-1))

        if( .not.(output.eq.TEMP_TO_ENT) ) then
	  xtemp = rtnewt_findtheta(ent_mismatch,tmin,tmax,RTACC,xrho,xye,xs)
	end if

        call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)

        xrho = xrho/rhounit
        xenr = xenr/eomunit
        xcs2 = xcs2/c2unit
        xprs = xprs/punit
	xent = xent/entunit
	!xent = xent*9.56565348d17/eomunit

        GammaP = log(xprs/kappa)/log(xrho)

	if(1.0+xenr+xprs/xrho.lt.hmin) then
	  hmin = 1.0+xenr+xprs/xrho
	end if

	  cs = 0.d0
	  if(xcs2.gt.0.d0) cs = sqrt(xcs2/(1.0+xenr+xprs/xrho))

        if(output.eq.MICRO) then
          write(6,"(1P10E15.6)") xrho,xxn,xxp,xxa,xxh	
	else if(output.eq.MUX) then
	  write(6,"(1P10E15.6)") xmu_e,xmu_p - xmu_n,xxn,xxp,xxh,xabar,xzbar
	else if(output.eq.EOSTABLE) then
          write(6,"(1P10E15.8)") xenr, GammaP, cs
	else if(output.eq.ENT_TO_TEMP) then
	  write(6,"(1P10E15.8)") xtemp
	else if(output.eq.TEMP_TO_ENT) then
	  write(6,"(1P10E15.8)") xent*sqrt(xrho)
	end if
      end do
    end do
  end do

!  write(*,*) "hmin=",hmin

end program driver


SUBROUTINE ent_mismatch(xtemp,fval,fderiv,xrho,xye,xentGoal)
  !USE nrtype
  IMPLICIT NONE
  REAL*8, INTENT(IN)  :: xrho,xye,xtemp,xentGoal
  REAL*8, INTENT(OUT) :: fval,fderiv
  real*8 :: xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe
  real*8 :: xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  integer keytemp,keyerr
  real*8 :: temp_p, temp_m, EPS, f_p, f_m, deriv_crit, tmin
  real*8 xunit,tunit,munit,rhounit,punit,eomunit,c2unit,entunit
  EPS = 1.d-5
  keytemp = 1
  keyerr  = 0
  xunit = 1.4766d5
  tunit = 4.9255d-6
  munit = 1.989d33
  rhounit = munit/xunit**3
  entunit = eomunit/9.56565348d17

  call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  fval =xent*sqrt(xrho/rhounit) - xentGoal

  temp_p = xtemp + EPS
  call nuc_eos_full(xrho,temp_p,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  f_p = xent*sqrt(xrho/rhounit) - xentGoal

  temp_m = xtemp - EPS
  call nuc_eos_full(xrho,temp_m,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr)
  f_m = xent*sqrt(xrho/rhounit) - xentGoal

  fderiv = (f_p - f_m)/(2.d0*EPS)
  

!  tmin = 0.1d0
!  deriv_crit = 1.01d0*fval/(xtemp-tmin)
!  if(abs(fderiv).lt.abs(deriv_crit)) then
!    fderiv = deriv_crit
!  end if


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
	  if(xent.lt.2.d0) then
	    rtnewt_findtheta = xl
	  else
	    rtnewt_findtheta = xh
	  end if
!	  write(*,*) "rtnewt_findtheta bracketing error",xent,fl,fh
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
