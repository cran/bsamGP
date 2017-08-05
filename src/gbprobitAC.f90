subroutine gbprobitAC(yobs,dobs,init_beta,beta_m0,beta_v0,nobs,npar, &
                      nburn,nthin,npost,ndisp,betaps,loglikeps,logpriorps)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nobs,npar,nburn,nthin,npost,ndisp
integer,intent(in) :: yobs(nobs)
real(8), intent(in) :: dobs(nobs,npar)
real(8), intent(in) :: init_beta(npar),beta_m0(npar),beta_v0(npar,npar)

! output arguments
real(8),intent(out) :: betaps(npost,npar),loglikeps(npost),logpriorps(npost)

! internal arguments
integer :: iobs
integer :: imcmc,isave,nmcmc
real(8)  :: stime,itime

real(8) :: dobst(npar,nobs),beta(npar),beta_v0i(npar,npar),beta_lnv0,loglik
real(8) :: Dbeta(nobs),beta_v0im0(npar),dtd(npar,npar),zobs(nobs)
real(8) :: beta_mn(npar),beta_vni(npar,npar),beta_vn(npar,npar)
real(8) :: cdfnorm,ltnormrnd,rtnormrnd

call cpu_time(stime)
call rndstart()

dobst=transpose(dobs)
dtd=matmul(dobst,dobs)

call inverse(beta_v0,npar,beta_v0i)
beta_v0im0=matmul(beta_v0i,beta_m0)
beta_lnv0=dlog(determinant(beta_v0,npar))
beta_vni=beta_v0i+dtd
call inverse(beta_vni,npar,beta_vn)

beta=init_beta
Dbeta=matmul(dobs,beta)

nmcmc=nburn+nthin*npost
isave=1
call dblepr('Burnin ...',-1,1.d0,0)
do imcmc=1,nmcmc
  if(imcmc.eq.nburn+1) call dblepr('Main iterations ...',-1,1.d0,0)

  call rchkusr() ! check interrupt

  do iobs=1,nobs
    if(yobs(iobs).eq.1) then
      zobs(iobs)=ltnormrnd(Dbeta(iobs),1.d0,0.d0)
    else
      zobs(iobs)=rtnormrnd(Dbeta(iobs),1.d0,0.d0)
    end if
  end do

  beta_mn=matmul(beta_vn,beta_v0im0+matmul(dobst,zobs))
  call mvnrnd(beta_mn,beta_vn,npar,beta)
  Dbeta=matmul(dobs,beta)

  if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
    betaps(isave,:)=beta
    logpriorps(isave)=GetLogPrior()
    loglik=0.d0
    do iobs=1,nobs
      if(yobs(iobs).eq.1) then
        loglik=loglik+cdfnorm(Dbeta(iobs),0.d0,1.d0,1,1)
      else
        loglik=loglik+cdfnorm(Dbeta(iobs),0.d0,1.d0,0,1)
      end if
    end do
    loglikeps(isave)=loglik

    if (mod(isave,ndisp).eq.0) then
      call cpu_time(itime)
      call sprint(isave,npost,itime-stime)
    end if

    isave=isave+1
  end if
end do
call rndend()

contains
!=========================================================================================

function GetLogPrior()
implicit none

!Output argument
real(8) :: GetLogPrior

GetLogPrior=-dot_product(beta-beta_m0,matmul(beta_v0i,beta-beta_m0))/2.d0- &
            dble(npar)*dlog(2.d0*PI)/2.d0-beta_lnv0/2.d0

return
end function GetLogPrior

end subroutine gbprobitAC
