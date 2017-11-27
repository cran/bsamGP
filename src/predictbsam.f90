subroutine predictbsam(xobs,xmin,xmax,nobs,nfun,nbasis,nint,fmodel,fpm,smcmc,&
                       thetag,alphag,psig,omegag,fxobsg)
use bsamTools
implicit none

!input arguments
integer,intent(in) :: nobs,nfun,nbasis,nint,fmodel(nfun),smcmc
real(8),intent(in) :: xobs(nobs,nfun),fpm(nfun),xmin(nfun),xmax(nfun)
real(8),intent(in) :: thetag(nbasis+1,nfun,smcmc),alphag(smcmc,nfun)
real(8),intent(in) :: psig(smcmc,nfun),omegag(smcmc,nfun)

!output arguments
real(8),intent(out) :: fxobsg(nobs,nfun,smcmc)

!internal arguments
integer :: imcmc,iloop,ifun,nr
integer :: quadfacts((nbasis+1)*(nbasis+2)/2,3),multfact((nbasis+1)*(nbasis+2)/2)
integer :: rightindex((nbasis+1)*(nbasis+2)/2),leftindex((nbasis+1)*(nbasis+2)/2)
integer :: tmpA(nbasis+1,1),tmpB(1,nbasis+1),tmpC(nbasis+1,nbasis+1)
integer :: tmpI(nbasis+1,nbasis+1),intsimpfacts(nint+1),xinxgrid(nobs,nfun)

real(8) :: xidelta(nobs,nfun),xdelta(nfun)
real(8) :: xrange(nfun),xgrid(nint+1,nfun)
real(8) :: xobs2(nobs),xgrid2(nint+1),xmid(nfun)

integer :: nfree,nfunconstraint,ifree,irest
real(8),allocatable :: phixobsfree(:,:,:),phixgridfree(:,:,:)
real(8),allocatable :: phixobs(:,:,:),phixgrid(:,:,:)

real(8) :: theta(nbasis+1,nfun)
real(8) :: psi(nfun),omega(nfun),alpha(nfun)
real(8) :: fxobs(nobs,nfun),fxgrid(nint+1,nfun)

! factors for integration
intsimpfacts(1)=1
intsimpfacts(nint+1)=1
do iloop=2,nint,2
  intsimpfacts(iloop)=4
end do
do iloop=3,nint,2
  intsimpfacts(iloop)=2
end do

tmpC=2
call Idiag(1,nbasis+1,tmpI)
call Ivech(tmpC-tmpI,nbasis+1,nbasis+1,multfact)
quadfacts(:,1)=multfact

tmpB=1
tmpA(:,1)=(/ (iloop,iloop=1,nbasis+1) /)
call Ikron(tmpB,1,nbasis+1,tmpA,nbasis+1,1,tmpC)
call Ivech(tmpC,nbasis+1,nbasis+1,leftindex)
quadfacts(:,2)=leftindex

tmpA=1
tmpB(1,:)=(/ (iloop,iloop=1,nbasis+1) /)
call Ikron(tmpA,nbasis+1,1,tmpB,1,nbasis+1,tmpC)
call Ivech(tmpC,nbasis+1,nbasis+1,rightindex)
quadfacts(:,3)=rightindex

! gridpoints for integration
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
do ifun=1,nfun
  if(fmodel(ifun).eq.5 .or. fmodel(ifun).eq.6 .or. fmodel(ifun).eq.7) then
    call intxgrid(xobs(:,ifun),xmin(ifun),xmax(ifun),xgrid(:,ifun),nobs,nint, &
                  xinxgrid(:,ifun),xidelta(:,ifun))
  end if
end do

! basis functions
nr=(nbasis+1)*(nbasis+2)/2
nfree=count(fmodel.eq.1)
nfunconstraint=count(fmodel.gt.1)
allocate(phixobsfree(nobs,nbasis,nfree),phixobs(nr,nobs,nfunconstraint))
allocate(phixgridfree(nint+1,nbasis,nfree),phixgrid(nr,nint+1,nfunconstraint))

ifree=1
irest=1
do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobsfree(:,:,ifree))
    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgridfree(:,:,ifree))
    ifree=ifree+1
  else
    if (fmodel(ifun).eq.2) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
                  IntConst2,IntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
                  IntConst2,IntCos,nbasis,nint+1,phixgrid(:,:,irest))
    else if (fmodel(ifun).eq.3) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))
    else if (fmodel(ifun).eq.4) then
      xobs2=xmin(ifun)+xmax(ifun)-xobs(:,ifun)
      xgrid2=xmin(ifun)+xmax(ifun)-xgrid(:,ifun)
      call GetPhi(xobs2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))
    else if (fmodel(ifun).eq.5) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
    else if (fmodel(ifun).eq.6) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
    else if (fmodel(ifun).eq.7) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
    end if
    irest=irest+1
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

  ifree=1
  irest=1
  do ifun=1,nfun
    if (fmodel(ifun).eq.1) then
      call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
                    nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))

      ifree=ifree+1
    else
      if (fmodel(ifun).eq.2) then
        call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest), &
                    quadfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.3) then
        call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                        xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                        nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.4) then
        call GetConcavef(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                         xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                         nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.5) then
        call GetSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),&
                   xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
                   xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
                   intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      else if (fmodel(ifun).eq.6) then
        call GetRotateSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun), &
                         xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
                         xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
                         xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
                         fxobs(:,ifun),fxgrid(:,ifun))
      else
        call GetUf(fpm(ifun),omega(ifun),psi(ifun),theta(:,ifun),xobs(:,ifun), &
                   xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
                   xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),quadfacts, &
                   intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
      end if

      irest=irest+1
    end if
  end do
  fxobsg(:,:,imcmc)=fxobs
end do
deallocate(phixobsfree,phixobs,phixgridfree,phixgrid)
end subroutine predictbsam
