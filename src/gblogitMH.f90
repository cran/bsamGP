subroutine gblogitMH(yint,xdata,init_beta,b,B0,nobs,nparx,&
                     nburn,nthin,nsave,ndisp,betas,loglikeps,logpriorps)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nobs,nparx,nburn,nthin,nsave,ndisp
integer,intent(in) :: yint(nobs)
real(8), intent(in) :: xdata(nobs,nparx),init_beta(nparx)
real(8), intent(in) :: b(nparx),B0(nparx,nparx)

! output arguments
real(8),intent(out) :: betas(nsave,nparx),loglikeps(nsave),logpriorps(nsave)

! internal arguments
integer :: imcmc,isave,nmcmc
real(8)  :: stime,itime

real(8) :: Xb(nobs),beta(nparx),iB0(nparx,nparx),yobs(nobs)
real(8) :: Ix(nparx,nparx),beta_vn(nparx,nparx)
real(8) :: beta_mAM(nparx,1),beta_vAM(nparx,nparx)

yobs=dble(yint)

! fixed parameters
call diag(1.d0,nparx,Ix)

! priors
call inverse(B0,nparx,iB0)

! initialization
beta=init_beta
Xb=matmul(xdata,beta)

! adaptation
beta_mAM(:,1)=beta

! Gibbs sampling
call cpu_time(stime)
call rndstart()
nmcmc=nburn+nthin*nsave
isave=1
do imcmc=1,nmcmc
  call rchkusr() ! check interrupt

  call update_beta()  ! update beta

  !save states
  if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
    ! save current state
    betas(isave,:)=beta

    logpriorps(isave)=mvnpdf(beta,b,B0,nparx,.true.)
    loglikeps(isave)=loglik_logit(yobs,Xb,nobs)

    if (mod(isave,ndisp).eq.0) then
      call cpu_time(itime)
      call sprint(isave,nsave,itime-stime)
    end if

    isave=isave+1
  end if
end do

call rndend()
!=========================================================================================

contains
!=========================================================================================


!-----------------------------------------------------------------------------------------
! *   Adaptive Metropolis (Roberts and Rosenthal, 2009, JCGS)
!-----------------------------------------------------------------------------------------
subroutine update_beta()
implicit none

!internal arguments
real(8) :: Beta0(nparx),Xb0(nobs),Beta1(nparx),Xb1(nobs)
real(8) :: loglik0,loglik1,beta_testp,beta_m1AM(nparx,1)
real(8) :: betat(1,nparx),rndunif

Beta0=beta
Xb0=Xb
loglik0=loglik_logit(yobs,Xb0,nobs)

if(imcmc.le.nparx*2) then
  beta_vn=(0.01d0)*Ix/dble(nparx)
else
  if(rndunif().le.0.05d0) then
    beta_vn=(0.01d0)*Ix/dble(nparx)
  else
    beta_vn=((2.38d0)**2.d0)*beta_vAM/dble(nparx)
  end if
end if
call mvnrnd(Beta0,beta_vn,nparx,Beta1)
Xb1=matmul(xdata,Beta1)
loglik1=loglik_logit(yobs,Xb1,nobs)

beta_testp=loglik1+mvnpdf(Beta1,b,B0,nparx,.true.)- &
loglik0-mvnpdf(Beta0,b,B0,nparx,.true.)
if(dlog(rndunif()).le.beta_testp) then
  beta=Beta1
  Xb=matmul(xdata,beta)
end if

! calculate empirical covariance matrix
beta_m1AM=beta_mAM
beta_mAM(:,1)=(dble(imcmc)*beta_m1AM(:,1)+beta)/dble(imcmc+1)
if(imcmc.eq.nparx*2) call covariance(betas(1:(imcmc+1),:),imcmc+1,nparx,beta_vAM)
if(imcmc.gt.nparx*2) then
  betat(1,:)=beta
  beta_vAM=(dble(imcmc-1)/dble(imcmc))*beta_vAM+ &
           (dble(imcmc)*matmul(beta_m1AM,transpose(beta_m1AM))- &
           dble(imcmc-1)*matmul(beta_mAM,transpose(beta_mAM))+ &
           matmul(transpose(betat),betat))/dble(imcmc)
end if

return
end subroutine update_beta


real(8) function loglik_logit(y,xb,nobs)
implicit none

!input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: y(nobs),xb(nobs)

loglik_logit=sum(y*xb-dlog(1.d0+dexp(xb)))

return
end function loglik_logit

end subroutine gblogitMH
