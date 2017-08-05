subroutine blreg(sigma2g,beta_mn,beta_vn,nparw,nmcmc,betag)
  use ToolsRfunf
  implicit none

  integer,intent(in) :: nparw,nmcmc
  real(8), intent(in) :: sigma2g(nmcmc),beta_mn(nparw),beta_vn(nparw,nparw)

  real(8),intent(out) :: betag(nmcmc,nparw)

  integer :: imcmc

  call rndstart()
  do imcmc=1,nmcmc
    call mvnrnd(beta_mn,sigma2g(imcmc)*beta_vn,nparw,betag(imcmc,:))
  end do
  call rndend()

  return
end subroutine blreg
