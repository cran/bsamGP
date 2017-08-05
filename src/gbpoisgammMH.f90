subroutine gbpoisgammMH(yint,xdata,init_beta,b,B0,kappa_m0,kappa_v0,nobs,nparx,&
                        nburn,nthin,nsave,ndisp,betas,kappas,loglikeps,logpriorps)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nobs,nparx,nburn,nthin,nsave,ndisp
integer,intent(in) :: yint(nobs)
real(8), intent(in) :: xdata(nobs,nparx),init_beta(nparx)
real(8), intent(in) :: b(nparx),B0(nparx,nparx),kappa_m0,kappa_v0

! output arguments
real(8),intent(out) :: betas(nsave,nparx),kappas(nsave)
real(8),intent(out) :: loglikeps(nsave),logpriorps(nsave)

! internal arguments
integer :: imcmc,isave,nmcmc
real(8)  :: stime,itime

integer :: nbatch,ibatch,naccept
real(8) :: yobs(nobs),Xb(nobs),beta(nparx),kappa,iB0(nparx,nparx),lambda(nobs)
real(8) :: Ix(nparx,nparx),beta_vn(nparx,nparx),r0,s0
real(8) :: beta_mAM(nparx,1),beta_vAM(nparx,nparx),lsd
real(8) :: dgamm

yobs=dble(yint)

call diag(1.d0,nparx,Ix)

call inverse(B0,nparx,iB0)
r0=2.d0*(2.d0+kappa_m0**2.d0/kappa_v0)
s0=kappa_m0*(r0-2.d0)

beta=init_beta
Xb=matmul(xdata,beta)
kappa=1.d0

beta_mAM(:,1)=beta
lsd=0.d0
ibatch=0
nbatch=50
naccept=0

call cpu_time(stime)
call rndstart()

nmcmc=nburn+nthin*nsave
isave=1
do imcmc=1,nmcmc
  call rchkusr()

  call update_lambda()
  call update_beta()
  call update_kappa()

  if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
    betas(isave,:)=beta
    kappas(isave)=kappa

    logpriorps(isave)=mvnpdf(beta,b,B0,nparx,.true.)+dgamm(kappa,r0,1.d0/s0,1)
    loglikeps(isave)=loglik_negbin(yobs,dexp(Xb),kappa,nobs)

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

subroutine update_lambda()
implicit none

!internal arguments
integer :: i
real(8) :: gamrnd

do i=1,nobs
  lambda(i)=gamrnd(yobs(i)+kappa,1.d0/(kappa/exp(Xb(i))+1.d0))
end do

return
end subroutine update_lambda


subroutine update_beta()
implicit none

!internal arguments
real(8) :: Beta0(nparx),Xb0(nobs),Beta1(nparx),Xb1(nobs)
real(8) :: beta_testp,beta_m1AM(nparx,1)
real(8) :: betat(1,nparx),rndunif

Beta0=beta
Xb0=Xb

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

beta_testp=loglikprior_beta(Beta1,lambda,Xb1,kappa,b,iB0,nobs,nparx)- &
loglikprior_beta(Beta0,lambda,Xb0,kappa,b,iB0,nobs,nparx)
if(dlog(rndunif()).le.beta_testp) then
  beta=Beta1
  Xb=matmul(xdata,beta)
end if

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


real(8) function loglikprior_beta(beta,lambda,Xb,kappa,b,iB0,nobs,nparx)
implicit none

! input arguments
integer,intent(in) :: nparx,nobs
real(8), intent(in) :: beta(nparx),lambda(nobs),Xb(nobs),kappa
real(8), intent(in) :: b(nparx),iB0(nparx,nparx)

! internal arguments
real(8) :: btb(1,1),mtb(1,1)

btb=dot_product(beta,matmul(iB0,beta))
mtb=dot_product(b,matmul(iB0,beta))

loglikprior_beta=-0.5d0*btb(1,1)-mtb(1,1)-kappa*(sum(Xb)+sum(lambda*dexp(-Xb)))

return
end function loglikprior_beta


subroutine update_kappa()
implicit none

!internal arguments
real(8) :: kappa0,kappa1,kappa_testp,normrnd,rndunif

if(mod(imcmc,nbatch).eq.0) then
  ibatch=ibatch+1
  if(dble(naccept)/dble(imcmc).lt.0.44d0) then
    lsd=lsd-min(0.01d0,1.d0/dsqrt(dble(ibatch)))
  else if(dble(naccept)/dble(imcmc).gt.0.44d0) then
    lsd=lsd+min(0.01d0,1.d0/dsqrt(dble(ibatch)))
  end if
end if

kappa0=kappa
kappa1=normrnd(kappa0,dexp(lsd))

kappa_testp=loglikprior_kappa(kappa1,lambda,Xb,nobs,r0,s0)- &
loglikprior_kappa(kappa0,lambda,Xb,nobs,r0,s0)

if(dlog(rndunif()).le.kappa_testp) then
  kappa=kappa1
  naccept=naccept+1
end if

return
end subroutine update_kappa


real(8) function loglikprior_kappa(kappa,lambda,Xb,nobs,r0,s0)
implicit none

! input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: kappa,lambda(nobs),Xb(nobs),r0,s0

! internal arguments
real(8) :: gammaln

loglikprior_kappa=(dble(nobs)*kappa+r0-1.d0)*dlog(kappa)-dble(nobs)*gammaln(kappa)- &
kappa*(sum(Xb))+(kappa-1.d0)*sum(dlog(lambda))- &
kappa*(s0+sum(lambda*dexp(-Xb)))

return
end function loglikprior_kappa


real(8) function loglik_negbin(y,mu,kappa,nobs)
implicit none

!input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: y(nobs),mu(nobs),kappa

!internal argument
integer :: i
real(8) :: loglik,const(nobs),gammaln

do i=1,nobs
  const(i)=gammaln(y(i)+kappa)-gammaln(y(i)+1.d0)
end do
loglik=-sum(y*dlog(kappa/mu+1.d0))-kappa*sum(dlog(1.d0+mu/kappa))- &
       dble(nobs)*gammaln(kappa)+sum(const)
loglik_negbin=loglik

return
end function loglik_negbin

end subroutine gbpoisgammMH
