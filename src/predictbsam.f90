subroutine predictbsam(xobs,xmin,xmax,nobs,nfun,nbasis,nint,fmodel,fpm,smcmc,&
                       thetag,alphag,psig,omegag,fxobsg)
use gbsamTools
implicit none

!input arguments
integer,intent(in) :: nobs,nfun,nbasis,nint,fmodel(nfun),smcmc
real(8),intent(in) :: xobs(nobs,nfun),fpm(nfun),xmin(nfun),xmax(nfun)
real(8),intent(in) :: thetag(nbasis+1,nfun,smcmc)
real(8),intent(in) :: alphag(smcmc,nfun),psig(smcmc,nfun),omegag(smcmc,nfun)

!output arguments
real(8),intent(out) :: fxobsg(nobs,nfun,smcmc)

!internal arguments
integer :: imcmc,iloop,ifun

integer :: xinxgrid(nobs,nfun),intsimpfacts(nint+1)
real(8) :: xobsint(nobs,nfun),xidelta(nobs,nfun)

real(8) :: xrange(nfun)
real(8) :: xgrid(nint+1,nfun),xmid(nfun),xdelta(nfun)

real(8) :: phixobs(nobs,nbasis+1,nfun)
real(8) :: phixgrid(nint+1,nbasis+1,nfun)

real(8) :: theta(nbasis+1,nfun)
real(8) :: psi(nfun),omega(nfun),alpha(nfun)
real(8) :: fxobs(nobs,nfun),fxgrid(nint+1,nfun)

intsimpfacts(1)=1
intsimpfacts(nint+1)=1
do iloop=2,nint,2
  intsimpfacts(iloop)=4
end do
do iloop=3,nint-1,2
  intsimpfacts(iloop)=2
end do

do ifun=1,nfun
  xmid(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
  xrange(ifun)=xmax(ifun)-xmin(ifun)
  xdelta(ifun)=(xrange(ifun))/dble(nint)
  xgrid(1,ifun)=xmin(ifun)
  do iloop=2,nint+1
    xgrid(iloop,ifun)=xgrid(iloop-1,ifun)+xdelta(ifun)
  end do
end do

xinxgrid=0
xidelta=0.d0
xobsint=xobs
do ifun=1,nfun
  if(fmodel(ifun).eq.4) then
    xobsint(:,ifun)=xmax(ifun)-xobs(:,ifun)+xmin(ifun)
  end if
  if(fmodel(ifun).gt.1) then
    call intxgrid(xobsint(:,ifun),xmin(ifun),xmax(ifun),xgrid(:,ifun),nobs,nint, &
                  xinxgrid(:,ifun),xidelta(:,ifun))
  end if
end do


do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobs(:,2:(nbasis+1),ifun))
    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgrid(:,2:(nbasis+1),ifun))
    phixobs(:,1,ifun)=0.d0
    phixgrid(:,1,ifun)=0.d0

    call GetFreef(theta(2:(nbasis+1),ifun),phixobs(:,2:(nbasis+1),ifun),phixgrid(:,2:(nbasis+1),ifun),&
                  nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
  else
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobs(:,2:(nbasis+1),ifun))
    call ConstFun(xobs(:,ifun),xrange(ifun),nobs,phixobs(:,1,ifun))

    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgrid(:,2:(nbasis+1),ifun))
    call ConstFun(xgrid(:,ifun),xrange(ifun),nint+1,phixgrid(:,1,ifun))
  end if
end do

do imcmc=1,smcmc
  call rchkusr() ! check interrupt

  fxobs=0.d0
  fxgrid=0.d0

  theta=thetag(:,:,imcmc)
  alpha=alphag(imcmc,:)
  psi=psig(imcmc,:)
  omega=omegag(imcmc,:)
  do ifun=1,nfun
    if (fmodel(ifun).eq.1) then
      call GetFreef(theta(2:(nbasis+1),ifun),phixobs(:,2:(nbasis+1),ifun),phixgrid(:,2:(nbasis+1),ifun),&
                    nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    else
      if (fmodel(ifun).eq.2) then
        call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun), &
                    xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun), &
                    intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.3) then
        call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                        xmid(ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                        xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                        intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.4) then
        call GetConcavef(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                         xmid(ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                         xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                         intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.5) then
        call GetSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),&
                   xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun), &
                   xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),&
                   intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.6) then
        call GetRotateSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun), &
                         xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                         xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),&
                         xmid(ifun),intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else
        call GetUf(fpm(ifun),omega(ifun),psi(ifun),theta(:,ifun),xobs(:,ifun), &
                   xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                   xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),&
                   intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      end if
    end if
  end do
  fxobsg(:,:,imcmc)=fxobs
end do

end subroutine predictbsam
