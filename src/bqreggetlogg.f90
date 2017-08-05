subroutine bqreggetlogg(betag,sigma2g,smcmc,npar,beta_mn,beta_covi,lndetbcov,sigma2_rn,sigma2_sn,logg)
  use ToolsRfunf
  implicit none

  integer,intent(in) :: npar,smcmc
  real(8), intent(in) :: betag(smcmc,npar),beta_mn(npar)
  real(8), intent(in) :: beta_covi(npar,npar),lndetbcov
  real(8), intent(in) :: sigma2g(smcmc),sigma2_rn,sigma2_sn

  real(8),intent(out) :: logg(smcmc)

  integer :: imcmc
  real(8) :: beta(npar),residbeta(npar),sigma2,logIlik,gammaln

  do imcmc=1,smcmc
    beta=betag(imcmc,:)
    sigma2=sigma2g(imcmc)

    logIlik=0.d0

    residbeta=beta-beta_mn
    logIlik=logIlik-dot_product(residbeta,matmul(beta_covi,residbeta))/2.d0- &
            dble(npar)*dlog(2.d0*PI)/2.d0-lndetbcov/2.d0

    logIlik=logIlik+(sigma2_rn/2.d0)*dlog(sigma2_sn/2.d0)-gammaln(sigma2_rn/2.d0)- &
            (sigma2_rn/2.d0+1.d0)*dlog(sigma2)-sigma2_sn/(2.d0*sigma2)

    logg(imcmc)=logIlik
  end do

  return
end subroutine bqreggetlogg
