! Surface PM from GOCART2G aerosol bins: PM = rho_air * sum_n( w_n * q_n ) at the lowest model level.
! w_n is the bin fraction below the PM size cut (and 1.38 for sulfate), as in GOCART2G_GridCompMod.F90.
! The operator is linear in the mixing ratios q_n. Density is taken from the trajectory and is not perturbed.
module PM_gocart_mod
  use aero_kinds_mod
  implicit none
  private
  public get_PM_gocart
  public get_PM_gocart_tl
  public get_PM_gocart_ad
  public find_lowest_level

contains

!> Index of the lowest model level of each profile: the level with the highest air pressure.
subroutine find_lowest_level(km, nobs, prs, ksfc)
  integer,              intent(in)  :: km
  integer,              intent(in)  :: nobs
  real(kind=kind_real), intent(in)  :: prs(km, nobs)
  integer,              intent(out) :: ksfc(nobs)
  integer :: i
  do i = 1, nobs
    if (prs(1, i) >= prs(km, i)) then
      ksfc(i) = 1
    else
      ksfc(i) = km
    end if
  end do
end subroutine find_lowest_level

!> Nonlinear: pm(iv,i) = rho(ksfc,i) * sum_n w(n) * qm(n,ksfc,i)
subroutine get_PM_gocart(km, nv, nobs, nq, w, ksfc, rho, qm, pm)
  integer,              intent(in)  :: km, nv, nobs, nq
  real(kind=kind_real), intent(in)  :: w(nq)
  integer,              intent(in)  :: ksfc(nobs)
  real(kind=kind_real), intent(in)  :: rho(km, nobs)
  real(kind=kind_real), intent(in)  :: qm(nq, km, nobs)
  real(kind=kind_real), intent(out) :: pm(nv, nobs)
  integer :: i, iq
  real(kind=kind_real) :: total
  do i = 1, nobs
    total = 0.0_kind_real
    do iq = 1, nq
      total = total + w(iq) * qm(iq, ksfc(i), i)
    end do
    pm(:, i) = rho(ksfc(i), i) * total
  end do
end subroutine get_PM_gocart

!> Tangent linear with trajectory density.
subroutine get_PM_gocart_tl(km, nv, nobs, nq, w, ksfc, rho, qm_tl, pm_tl)
  integer,              intent(in)  :: km, nv, nobs, nq
  real(kind=kind_real), intent(in)  :: w(nq)
  integer,              intent(in)  :: ksfc(nobs)
  real(kind=kind_real), intent(in)  :: rho(km, nobs)
  real(kind=kind_real), intent(in)  :: qm_tl(nq, km, nobs)
  real(kind=kind_real), intent(out) :: pm_tl(nv, nobs)
  call get_PM_gocart(km, nv, nobs, nq, w, ksfc, rho, qm_tl, pm_tl)
end subroutine get_PM_gocart_tl

!> Adjoint: qm_ad(n,ksfc,i) += rho * w(n) * sum_iv pm_ad(iv,i); other levels untouched.
subroutine get_PM_gocart_ad(km, nv, nobs, nq, w, ksfc, rho, pm_ad, qm_ad)
  integer,              intent(in)    :: km, nv, nobs, nq
  real(kind=kind_real), intent(in)    :: w(nq)
  integer,              intent(in)    :: ksfc(nobs)
  real(kind=kind_real), intent(in)    :: rho(km, nobs)
  real(kind=kind_real), intent(in)    :: pm_ad(nv, nobs)
  real(kind=kind_real), intent(inout) :: qm_ad(nq, km, nobs)
  integer :: i, iq
  real(kind=kind_real) :: pm_sum
  do i = nobs, 1, -1
    pm_sum = sum(pm_ad(:, i))
    do iq = nq, 1, -1
      qm_ad(iq, ksfc(i), i) = qm_ad(iq, ksfc(i), i) + rho(ksfc(i), i) * w(iq) * pm_sum
    end do
  end do
end subroutine get_PM_gocart_ad

end module PM_gocart_mod
