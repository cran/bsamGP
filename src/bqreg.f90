subroutine bqreg(yobs,dobs,p,init_beta,beta_m0,beta_v0,init_sigsq,sigma2_m0,sigma2_v0, &
                 nobs,npar,nburn,nthin,npost,betaps,sigmasqps,loglikeps,logpriorps)
  use ToolsRfunf
  implicit none

  integer,intent(in) :: nobs,npar,nburn,nthin,npost
  real(8), intent(in) :: yobs(nobs),dobs(nobs,npar),p
  real(8), intent(in) :: init_beta(npar),beta_m0(npar),beta_v0(npar,npar)
  real(8), intent(in) :: init_sigsq,sigma2_m0,sigma2_v0

  real(8),intent(out) :: betaps(npost,npar),sigmasqps(npost)
  real(8),intent(out) :: loglikeps(npost),logpriorps(npost)

  integer :: imcmc,isave,nmcmc

  real(8) :: beta(npar),nu(nobs),sigmasq,dobst(npar,nobs),beta_lnv0,beta_v0i(npar,npar)
  real(8) :: Dbeta(nobs),beta_v0im0(npar),ehat(nobs),eta1,eta1sq,eta2sq
  real(8) :: sigma2_r0,sigma2_s0,sigma2_rn,nu_v(nobs,nobs),nu_vi(nobs,nobs)

  dobst=transpose(dobs)

  eta1=(1.d0-2.d0*p)/(p*(1.d0-p))
  eta1sq=eta1**2.d0
  eta2sq=2.d0/(p*(1.d0-p))

  call inverse(beta_v0,npar,beta_v0i)
  beta_v0im0=matmul(beta_v0i,beta_m0)
  beta_lnv0=dlog(determinant(beta_v0,npar))
  sigma2_r0=2.d0*(2.d0+sigma2_m0**2.d0/sigma2_v0)
  sigma2_s0=sigma2_m0*(sigma2_r0-2.d0)
  sigma2_rn=sigma2_r0+3.d0*dble(nobs)+dble(npar)

  beta=init_beta
  Dbeta=matmul(dobs,beta)
  sigmasq=init_sigsq
  nu=sigmasq
  call diag(sigmasq,nobs,nu_v)
  call diag(sigmasq,nobs,nu_vi)

  call rndstart()
  nmcmc=nburn+nthin*npost
  isave=1
  do imcmc=1,nmcmc
    call rchkusr()

    call draw_beta()
    call draw_nu()
    call draw_sigma()

    if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
      betaps(isave,:)=beta
      sigmasqps(isave)=sigmasq
      logpriorps(isave)=GetLogPrior()
      ehat=yobs-Dbeta
      loglikeps(isave)=-sum(dabs(ehat)+(2.d0*p-1.d0)*ehat)/(2.d0*sigmasq)+ &
                       dble(nobs)*(dlog(p)+dlog(1.d0-p)-dlog(sigmasq))
      isave=isave+1
    end if
  end do
  call rndend()

contains

subroutine draw_beta()
  implicit none

  real(8) :: resid(nobs),vni(npar,npar),vn(npar,npar),bn(npar)

  resid=yobs-eta1*nu
  vni=matmul(dobst,matmul(nu_vi,dobs))/eta2sq+beta_v0i
  call inverse(vni,npar,vn)
  bn=beta_v0im0+matmul(dobst,matmul(nu_vi,resid))/eta2sq
  bn=matmul(vn,bn)

  call mvnrnd(bn,sigmasq*vn,npar,beta)
  Dbeta=matmul(dobs,beta)

  return
end subroutine draw_beta


subroutine draw_nu()
  implicit none

  integer :: iobs
  real(8) :: resid(nobs),lambda1(nobs),lambda2,invgaussrnd

  resid=yobs-Dbeta
  lambda2=2.d0/sigmasq+eta1sq/(sigmasq*eta2sq)
  lambda1=dsqrt((eta1sq+2.d0*eta2sq)/(resid**2.d0))
  do iobs=1,nobs
    nu(iobs)=1.d0/invgaussrnd(lambda1(iobs),lambda2)
  end do
  call diagvec(nu,nobs,nu_v)
  call diagvec(1.d0/nu,nobs,nu_vi)

  return
end subroutine draw_nu


subroutine draw_sigma()
  implicit none

  ! Internal arguments
  real(8) :: resid(nobs),residb(npar),sigma2_sn,gamrnd

  sigma2_sn=0.d0

  resid=yobs-Dbeta-eta1*nu
  sigma2_sn=sigma2_sn+sigma2_s0+sum(resid**2.d0/nu)/eta2sq

  sigma2_sn=sigma2_sn+2.d0*sum(nu)

  residb=beta-beta_m0
  sigma2_sn=sigma2_sn+dot_product(residb,matmul(beta_v0i,residb))

  sigmasq=1.d0/gamrnd(sigma2_rn/2.d0,2.d0/sigma2_sn)

  return
end subroutine draw_sigma


function GetLogPrior()
  implicit none

  real(8) :: GetLogPrior

  real(8) :: residb(npar),lnpriorf
  real(8) :: gammaln

  lnpriorf=0.d0

  residb=beta-beta_m0
  lnpriorf=lnpriorf-dot_product(residb,matmul(beta_v0i,residb))/(2.d0*sigmasq)- &
  dble(npar)*dlog(2.d0*PI*sigmasq)/2.d0-beta_lnv0/2.d0

  lnpriorf=lnpriorf+(sigma2_r0/2.d0)*dlog(sigma2_s0/2.d0)-gammaln(sigma2_r0/2.d0)- &
           (sigma2_r0/2.d0+1.d0)*dlog(sigmasq)-sigma2_s0/(2.d0*sigmasq)
  GetLogPrior=lnpriorf

  return
end function GetLogPrior

end subroutine bqreg
