subroutine gbglmgetlogg(betag,smcmc,npar,beta_mn,beta_covi,lndetbcov,logg)
use ToolsRfunf
implicit none

!input arguments
integer,intent(in) :: npar,smcmc
real(8), intent(in) :: betag(smcmc,npar),beta_mn(npar)
real(8), intent(in) :: beta_covi(npar,npar),lndetbcov

!output argument
real(8),intent(out) :: logg(smcmc)

!internal arguments
integer :: imcmc
real(8) :: beta(npar),residbeta(npar),logIlik

do imcmc=1,smcmc
beta=betag(imcmc,:)

logIlik=0.d0

residbeta=beta-beta_mn
logIlik=logIlik-dot_product(residbeta,matmul(beta_covi,residbeta))/2.d0- &
        dble(npar)*dlog(2.d0*PI)/2.d0-lndetbcov/2.d0

logg(imcmc)=logIlik
end do

return
end subroutine gbglmgetlogg
