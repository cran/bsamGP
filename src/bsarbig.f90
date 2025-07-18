subroutine bsarbig(verbose,ydata,xdata,nobs,nxgrid,nbasis,&
                   w0,tau2_m0,tau2_v0,sigma2_r0,sigma2_s0,&
                   nblow,smcmc,nskip,ndisp,fhatg,thetag,sigmag,smoothg)
use ToolsRfunf
implicit none

! input arguments
integer,intent(in) :: nbasis,nxgrid,nobs,verbose
integer,intent(in) :: nblow,smcmc,nskip,ndisp
real(8),intent(in) :: ydata(nobs),xdata(nobs),w0
real(8),intent(in) :: tau2_m0,tau2_v0,sigma2_r0,sigma2_s0

! output arguments
real(8),intent(out) :: fhatg(smcmc,nxgrid),thetag(smcmc,nbasis)
real(8),intent(out) :: sigmag(smcmc),smoothg(smcmc,2)

! internal arguments
real(8) :: stime,itime
integer :: imcmc,isave,nmcmc

integer :: i,j,k,nmiss
integer :: ndata(nxgrid),ndata0(nxgrid),id_xmiss(nxgrid)
real(8) :: kall(nbasis),phi(nxgrid,nbasis),nmat(nxgrid,nxgrid)
real(8) :: ptp(nbasis,nbasis),pty(nbasis,1),theta(nbasis,1),theta2(nbasis)
real(8) :: sigma,sigma2,tau,tau2,gampar,egamkall(nbasis),yhat(nxgrid,1)
real(8) :: wk,xdelta,xgrid(nxgrid),a(nxgrid),ybar(nxgrid,1)
real(8) :: vni(nbasis,nbasis),vn(nbasis,nbasis)
real(8) :: thetan(nbasis,1),resid(nxgrid),sse

real(8) :: tau2_r0,tau2_s0
real(8) :: sigma2_rn,sigma2_sn
real(8) :: rndnorm,gamrnd

call cpu_time(stime)
call rndstart()

kall=(/ (dble(k),k=0,nbasis-1) /)
ndata=0

! compute grid on 0 to 1
xdelta=1.d0/dble(nxgrid)
do i=1,nxgrid
  xgrid(i)=xdelta/2.d0+dble(i-1)*xdelta
end do

! Get summary statistics on intervals aroung grid points
a=xgrid+xdelta/2.d0
ndata(1)=count(xdata.le.a(1))
if (ndata(1).gt.0) then
  ybar(1,1)=sum(ydata*transfer(xdata.le.a(1),(/(1,i=1,nobs)/)))/dble(ndata(1))
end if
ndata(nxgrid)=count(xdata.gt.a(nxgrid-1))
if (ndata(nxgrid).gt.0) then
  ybar(nxgrid,1)=sum(ydata*transfer(xdata.gt.a(nxgrid-1),(/(1,i=1,nobs)/)))/dble(ndata(nxgrid))
end if
do j=2,nxgrid-1
  ndata(j)=count(xdata.gt.a(j-1) .and. xdata.le.a(j))
  if (ndata(j).gt.0) then
    ybar(j,1)=sum(ydata*transfer(xdata.gt.a(j-1) .and. xdata.le.a(j),(/(1,i=1,nobs)/)))/dble(ndata(j))
  end if
end do

! Check for intervals with 0 observation
call which(ndata.eq.0,nxgrid,id_xmiss,nmiss)
if (nmiss.gt.0) then
  ndata0=ndata
  ndata=ndata+transfer(ndata.eq.0,(/(1,i=1,nxgrid)/))
end if
call diagvec(dble(ndata),nxgrid,nmat)

! Compute basis function for density on xgrid
phi=1.d0
do j=2,nbasis
  do i=1,nxgrid
    phi(i,j)=dsqrt(2.d0)*cos(PI*xgrid(i)*kall(j))
  end do
end do
ptp=matmul(matmul(transpose(phi),nmat),phi)
pty=matmul(matmul(transpose(phi),nmat),ybar)

! Initialize priors for Fourier series
tau2_r0=2.d0*(2.d0+tau2_m0**2.d0/tau2_v0)
tau2_s0=tau2_m0*(tau2_r0-2.d0)
wk=sum(kall)/2.d0-w0

! Prior parameters
sigma2_rn=sigma2_r0+dble(nxgrid)+dble(nbasis)+dble(nmiss)

! Initialize parameters
theta=0.d0
sigma=1.d0
sigma2=sigma**2.d0
tau=1.d0
tau2=tau**2.d0
gampar=1.d0
egamkall=dexp(-gampar*kall)
yhat=matmul(phi,theta)

isave=1
nmcmc=nblow+nskip*smcmc
do imcmc=1,nmcmc
  if(imcmc.eq.1 .and. verbose.eq.1) then
    call biprint()
  end if
  if(imcmc.eq.nblow+1 .and. verbose.eq.1) then
    call miprint()
  end if

  call rchkusr()  ! user interrupt
  if(nmiss.gt.0) then
    do i=1,nmiss
      ybar(id_xmiss(i),1)=yhat(id_xmiss(i),1)+sigma*rndnorm()
    end do
    pty=matmul(matmul(transpose(phi),nmat),ybar)
  end if

  ! generate theta
  vni=ptp/sigma2
  do k=1,nbasis
    vni(k,k)=vni(k,k)+1.d0/(tau2*egamkall(k))
  end do
  call inverse(vni,nbasis,vn)
  thetan=matmul(vn,pty/sigma2)
  call mvnrnd(thetan(:,1),vn,nbasis,theta(:,1))
  theta2=theta(:,1)**2.d0

  ! generate smoothing parameters
  call getsmooth(theta2,nbasis,wk,w0,sigma2,tau2_r0,tau2_s0,tau2,gampar,egamkall)
  tau=dsqrt(tau2)

  ! generate sigma2
  yhat=matmul(phi,theta)
  resid=ybar(:,1)-yhat(:,1)
  sse=sum((resid**2.d0)*dble(ndata))+sum(theta2/(tau2*egamkall))
  sigma2_sn=sigma2_s0+sse
  sigma2=sigma2_sn/(2.d0*gamrnd(sigma2_rn/2.d0,1d0))
  sigma=dsqrt(sigma2)

  if(imcmc.gt.nblow .and. mod(imcmc,nskip).eq.0) then
    ! save posteriors
    fhatg(isave,:)=yhat(:,1)
    thetag(isave,:)=theta(:,1)
    sigmag(isave)=sigma
    smoothg(isave,1)=tau
    smoothg(isave,2)=gampar

    if (verbose.eq.1) then
      if (mod(isave,ndisp).eq.0) then
        call cpu_time(itime)
        call sprint(isave,smcmc,itime-stime)
      end if
    end if
    isave=isave+1
  end if
end do

call rndend()

!===============================================================================
contains


subroutine getsmooth(theta2,nbasis,wk,w0,sigma2,tau2_r0,tau2_s0,&
                     tau2,gampar,egamkall)
implicit none

! input arguments
integer,intent(in) :: nbasis
real(8),intent(in) :: theta2(nbasis),wk,w0,sigma2
real(8),intent(in) :: tau2_r0,tau2_s0

! output arguments
real(8),intent(out) :: tau2
real(8),intent(inout) :: gampar,egamkall(nbasis)

! internal arguments
integer :: k,iloop,ntheta0,z(nbasis)
real(8) :: kall(nbasis),tau2_rn,tau2_sn,ck(nbasis)
real(8) :: u1(nbasis),ckz(nbasis),u2,bnz(nbasis),bmin,rndunif
real(8),allocatable :: bnz1(:)

kall=(/ (dble(k),k=1,nbasis) /)
tau2_rn=tau2_r0+dble(nbasis)
tau2_sn=tau2_s0+sum(theta2/egamkall)/sigma2
tau2=tau2_sn/(2.d0*gamrnd(tau2_rn/2.d0,1.d0))

ck=theta2/(2.d0*tau2*sigma2)
ntheta0=count(ck.eq.0.d0)
if (ntheta0.eq.nbasis) then
  gampar = 1.d0/w0
else
  z=merge(1,0,ck.eq.0.d0)
  ckz=dble(z)*(1.d0-ck)+ck
  u1=(/ (rndunif(),k=1,nbasis) /)
  bnz=gampar+(dlog(ckz-dlog(u1)*egamkall)-dlog(ckz))/kall
  if (ntheta0.gt.0) then
    allocate(bnz1(count(z.ne.1)))
    iloop=1
    do k=1,nbasis
      if (z(k).ne.1) then
        bnz1(iloop)=bnz(k)
        iloop=iloop+1
      end if
    end do
    bmin=minval(bnz1)
    deallocate(bnz1)
  else
    bmin=minval(bnz)
  end if
  u2=rndunif()
  gampar=bmin+dlog(u2+(1.d0-u2)*dexp(-wk*bmin))/wk
end if
egamkall=dexp(-gampar*kall)

return
end subroutine getsmooth

end subroutine bsarbig
