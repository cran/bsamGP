subroutine gbpoisMH(verbose,yint,xdata,init_beta,b,B0,mub,Sb,nobs,nparx,&
                    nburn,nthin,nsave,ndisp,betas,loglikeps,logpriorps)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nobs,nparx,nburn,nthin,nsave,ndisp
integer,intent(in) :: yint(nobs),verbose
real(8), intent(in) :: xdata(nobs,nparx),init_beta(nparx)
real(8), intent(in) :: b(nparx),B0(nparx,nparx)
real(8), intent(in) :: mub(nparx),Sb(nparx,nparx)

! output arguments
real(8),intent(out) :: betas(nsave,nparx),loglikeps(nsave),logpriorps(nsave)

! internal arguments
integer :: iobs,imcmc,isave,nmcmc
real(8)  :: stime,itime

real(8) :: Xb(nobs),yobs(nobs)
real(8) :: beta(nparx),iB0(nparx,nparx),const(nobs)
real(8) :: gammaln

yobs=dble(yint)

call inverse(B0,nparx,iB0)

beta=init_beta
Xb=matmul(xdata,beta)

call cpu_time(stime)
call rndstart()

nmcmc=nburn+nthin*nsave
isave=1
if (verbose.eq.1) then
  call dblepr('Burnin ...',-1,1.d0,0)
end if
do imcmc=1,nmcmc
  call rchkusr() ! check interrupt
  if(imcmc.eq.nburn+1 .and. verbose.eq.1) call dblepr('Main iterations ...',-1,1.d0,0)

  call update_beta()

  if(imcmc .gt. nburn .and. mod(imcmc,nthin) .eq. 0) then
    betas(isave,:)=beta

    logpriorps(isave)=mvnpdf(beta,b,B0,nparx,.true.)
    do iobs=1,nobs
      const(iobs)=gammaln(yobs(iobs)+1.d0)
    end do
    loglikeps(isave)=-sum(dexp(Xb))+sum(yobs*Xb)-sum(const)

    if (verbose.eq.1) then
      if (mod(isave,ndisp).eq.0) then
        call cpu_time(itime)
        call sprint(isave,nsave,itime-stime)
      end if
    end if

    isave=isave+1
  end if
end do

call rndend()
!=========================================================================================

contains
!=========================================================================================

subroutine update_beta()
implicit none

!internal arguments
real(8) :: Beta0(nparx),Xb0(nobs),Beta1(nparx),Xb1(nobs)
real(8) :: lacc,rndunif

Beta0=Beta
Xb0=Xb

call mvnrnd(mub,Sb,nparx,Beta1)
Xb1=matmul(xdata,Beta1)

lacc=accept_ftn(yobs,Beta0,Xb0,Beta1,Xb1,b,iB0,mub,Sb,nparx,nobs)
if(dlog(rndunif()).le.lacc) then
  Beta=Beta1
  Xb=matmul(xdata,Beta)
end if

return
end subroutine update_beta

function accept_ftn(y,Beta0,Xb0,Beta1,Xb1,b,iB0,mub,Sb,nparx,nobs)
implicit none

!input arguments
integer,intent(in) :: nparx,nobs
real(8), intent(in) :: y(nobs),Beta0(nparx),Xb0(nobs),Beta1(nparx),Xb1(nobs)
real(8), intent(in) :: b(nparx),iB0(nparx,nparx),mub(nparx),Sb(nparx,nparx)

!output argument
real(8) :: accept_ftn

!internal arguments
real(8) :: lprior0,llike0,lprop0,lprior1,llike1,lprop1, residb(nparx)

residb=Beta0-b
lprior0=-dot_product(residb,matmul(iB0,residb))/2.d0
llike0=sum(-dexp(Xb0))+sum(y*Xb0)
lprop0=mvnpdf(Beta0,mub,Sb,nparx,.true.)

residb=Beta1-b
lprior1=-dot_product(residb,matmul(iB0,residb))/2.d0
llike1=sum(-dexp(Xb1))+sum(y*Xb1)
lprop1=mvnpdf(Beta1,mub,Sb,nparx,.true.)

accept_ftn=llike1+lprior1-lprop1-llike0-lprior0+lprop0

return
end function accept_ftn

end subroutine gbpoisMH
