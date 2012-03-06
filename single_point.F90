program driver

  use eosmodule
  implicit none

  real*8 xrho,xye,xtemp
  real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
  real*8 xdpderho,xdpdrhoe
  integer keytemp,keyerr

  ! for full eos call:
  real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  real*8 xxa,xxh,xxn,xxp

  real*8 xunit,tunit,munit,rhounit,punit,eomunit,c2unit
  real*8 BARYMASS

  keytemp = 1
  keyerr  = 0

  xunit = 1.4766d5
  tunit = 4.9255d-6
  munit = 1.989d33
  BARYMASS = 1.659d-24
  rhounit = munit/xunit**3
  punit = munit/(xunit*tunit**2)
  c2unit = (xunit/tunit)**2
  eomunit = c2unit

!  xrho  = 1.25957d-6*rhounit
!  xye   = 0.146549d0
!  xtemp = 5.794424d0

!  xrho   = 1.469967d-6*rhounit
!  xye    = 0.166106
!  xtemp  = 15.33564

!   xrho   = 6.8d-7*rhounit
!   xtemp  = 6.28d0
!   xye    = 0.114d0

    xrho   = 1.393d-6*rhounit
    xtemp  = 1.6629d0
    xye    = 7.248d-2

!   write(*,*) "rhounit=",rhounit

    call readtable("LS220_450r_270t_50y_062211.h5")

        call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
          xmuhat,keytemp,keyerr)

!   write(*,*) "eta_e=",xmu_e/xtemp
    write(*,*) "ent=",xent*9.56565348d17/eomunit

end program driver

