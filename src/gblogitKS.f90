subroutine gblogitKS(verbose,yobs,dobs,init_beta,beta_m0,beta_v0,nobs,npar, &
                     nburn,nthin,npost,ndisp,betaps,loglikeps,logpriorps)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nobs,npar,nburn,nthin,npost,ndisp
integer,intent(in) :: yobs(nobs),verbose
real(8), intent(in) :: dobs(nobs,npar)
real(8), intent(in) :: init_beta(npar),beta_m0(npar),beta_v0(npar,npar)

! output arguments
real(8),intent(out) :: betaps(npost,npar),loglikeps(npost),logpriorps(npost)

! internal arguments
integer :: iobs,ok
integer :: imcmc,isave,nmcmc
real(8)  :: stime,itime

real(8) :: dobst(npar,nobs),beta(npar),beta_v0i(npar,npar),beta_lnv0
real(8) :: Dbeta(nobs),beta_v0im0(npar),zobs(nobs),loglik,w(nobs,nobs)
real(8) :: beta_mn(npar),beta_vni(npar,npar),beta_vn(npar,npar)
real(8) :: r2(nobs),lamp,uacc,lam(nobs),dtwd(npar,npar),dtwz(npar)

real(8) :: rndunif,ltlogisrnd,rtlogisrnd

call cpu_time(stime)
call rndstart()

dobst=transpose(dobs)

call inverse(beta_v0,npar,beta_v0i)
beta_v0im0=matmul(beta_v0i,beta_m0)
beta_lnv0=dlog(determinant(beta_v0,npar))
lam=1.d0

beta=init_beta
Dbeta=matmul(dobs,beta)

nmcmc=nburn+nthin*npost
isave=1
if (verbose.eq.1) then
  call dblepr('Burnin ...',-1,1.d0,0)
end if
do imcmc=1,nmcmc
  if(imcmc.eq.nburn+1 .and. verbose.eq.1) call dblepr('Main iterations ...',-1,1.d0,0)

  do iobs=1,nobs
    call rchkusr() ! check interrupt

    if(yobs(iobs).eq.1) then
      zobs(iobs)=ltlogisrnd(Dbeta(iobs),1.d0,0.d0)
    else
      zobs(iobs)=rtlogisrnd(Dbeta(iobs),1.d0,0.d0)
    end if
  end do

  r2=(zobs-Dbeta)**2.d0
  do iobs=1,nobs
    call rchkusr() ! check interrupt

    lamp=gigrnd(0.5d0,1.d0,r2(iobs))
    uacc=rndunif()
    if(lamp.gt.(4.d0/3.d0)) then
      ok=rightmost_interval(uacc,lamp)
    else
      ok=leftmost_interval(uacc,lamp)
    end if
    if (ok.eq.1) lam(iobs)=lamp
  end do

  call diagvec(1.d0/lam,nobs,w)
  dtwd=matmul(dobst,matmul(w,dobs))
  beta_vni=beta_v0i+dtwd
  call inverse(beta_vni,npar,beta_vn)
  dtwz=matmul(dobst,matmul(w,zobs))
  beta_mn=matmul(beta_vn,beta_v0im0+dtwz)
  call mvnrnd(beta_mn,beta_vn,npar,beta)
  Dbeta=matmul(dobs,beta)

  if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
    betaps(isave,:)=beta
    logpriorps(isave)=GetLogPrior()
    loglik=0.d0
    do iobs=1,nobs
      loglik=loglik+dble(yobs(iobs))*Dbeta(iobs)-dlog(1.d0+dexp(Dbeta(iobs)))
    end do
    loglikeps(isave)=loglik

    if (verbose.eq.1) then
      if (mod(isave,ndisp).eq.0) then
        call cpu_time(itime)
        call sprint(isave,npost,itime-stime)
      end if
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

end subroutine gblogitKS
