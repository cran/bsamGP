module bsamTools
use ToolsRfunf
implicit none

contains

function LogfiG(x,r0,s0)
  implicit none

  !input arguments
  real(8),intent(in) :: x,r0,s0

  !output argument
  real(8) :: LogfiG

  !internal arguments
  real(8) :: gammaln

  LogfiG=(r0/2.d0)*dlog(s0/2.d0)-gammaln(r0/2.d0)- &
         (r0/2.d0+1.d0)*dlog(x)-s0/(2.d0*x)

  return
end function LogfiG


subroutine intxgrid(xobs,xmin,xmax,xgrid,nobs,nint,xinxgrid,xidelta)
  implicit none

  ! input arguments
  integer,intent(in) :: nobs,nint
  real(8), intent(in) :: xobs(nobs),xmin,xmax,xgrid(nint+1)

  ! output arguments
  integer,intent(out) :: xinxgrid(nobs)
  real(8), intent(out) :: xidelta(nobs)

  ! internal arguments
  integer :: iobs,iint,s(nint+1),szbot
  real(8) :: xi

  s=(/ (iint,iint=1,nint+1) /)

  do iobs=1,nobs
    xi=xobs(iobs)
    if(xi.eq.xmin) then
      xinxgrid(iobs)=1
    else if(xi.eq.xmax) then
      xinxgrid(iobs)=nint+1
    else
      szbot=0
      do iint=1,nint+1
        if(xgrid(iint).le.xi) szbot=iint
      end do
      xinxgrid(iobs)=szbot
      if (xi.gt.xgrid(szbot)) then
        xidelta(iobs)=xi-xgrid(szbot)
      end if
    end if
  end do

  return
end subroutine intxgrid


subroutine QuadMult(x,qvech,quadfacts,n,nr,nc,quadvec)
  implicit none

  !Input arguments
  integer,intent(in) :: n,nr,nc
  integer,intent(in) :: quadfacts(nr,3)
  real(8), intent(in) :: x(n),qvech(nr,nc)

  !Output argument
  real(8),intent(out) :: quadvec(nc)

  !Internal argument
  integer :: i
  real(8) :: tmp(nr,nc)

  do i=1,nc
    tmp(:,i)=quadfacts(:,1)*x(quadfacts(:,2))*qvech(:,i)*x(quadfacts(:,3))
  end do
  quadvec=sum(tmp,1)

  return
end subroutine QuadMult


subroutine GetFreef(theta,phixobs,phixgrid,nbasis,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nobs,ngrid
  real(8), intent(in) :: theta(nbasis),phixobs(nobs,nbasis),phixgrid(ngrid,nbasis)

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  fxobs=matmul(phixobs,theta)
  fxgrid=matmul(phixgrid,theta)

  return
end subroutine GetFreef


subroutine GetUpf(fpm,theta,phixobs,phixgrid,quadfacts,nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: quadfacts(nr,3)
  real(8), intent(in) :: fpm,theta(nbasis+1),phixobs(nr,nobs),phixgrid(nr,ngrid)

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,fxobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,fxgrid)

  fxgrid = fpm*fxgrid
  fxobs = fpm*fxobs

  return
end subroutine GetUpf


subroutine GetConvexf(fpm,alpha,theta,xobs,xgrid,xmid,phixobs,phixgrid,quadfacts, &
                      nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: quadfacts(nr,3)
  real(8), intent(in) :: fpm,alpha,theta(nbasis+1),xobs(nobs),xgrid(ngrid),xmid
  real(8), intent(in) :: phixobs(nr,nobs),phixgrid(nr,ngrid)

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,fxobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,fxgrid)

  fxgrid=fpm*fxgrid+alpha*(xgrid-xmid)
  fxobs=fpm*fxobs+alpha*(xobs-xmid)

  return
end subroutine GetConvexf


subroutine GetConcavef(fpm,alpha,theta,xobs,xgrid,xmid,phixobs,phixgrid,quadfacts, &
                       nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: quadfacts(nr,3)
  real(8), intent(in) :: fpm,alpha,theta(nbasis+1),xobs(nobs),xgrid(ngrid),xmid
  real(8), intent(in) :: phixobs(nr,nobs),phixgrid(nr,ngrid)

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,fxobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,fxgrid)

  fxgrid=-fxgrid
  fxobs=-fxobs

  fxgrid=fpm*fxgrid+alpha*(xgrid-xmid)
  fxobs=fpm*fxobs+alpha*(xobs-xmid)

  return
end subroutine GetConcavef


subroutine GetSf(fpm,omega,psi,alpha,theta,xobs,xgrid,phixobs,phixgrid,xdelta,xinxgrid, &
                 xidelta,xrange,xmid,quadfacts,intsimpfacts,nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: xinxgrid(nobs),quadfacts(nr,3),intsimpfacts(ngrid)
  real(8), intent(in) :: fpm,omega,psi,alpha,theta(nbasis+1),xobs(nobs),xgrid(ngrid)
  real(8), intent(in) :: phixobs(nr,nobs),phixgrid(nr,ngrid),xdelta,xidelta(nobs)
  real(8), intent(in) :: xrange,xmid

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  !internal arguments
  real(8) :: z2xobs(nobs),z2xgrid(ngrid),hxgrid(ngrid),hxobs(nobs)
  real(8) :: f2xgrid(ngrid),f2xobs(nobs),f1xgrid(ngrid),f1xobs(nobs),f1min,fint

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,z2xobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,z2xgrid)

  call SquishDown(xgrid,psi,omega,ngrid,hxgrid)
  call SquishDown(xobs,psi,omega,nobs,hxobs)

  f2xgrid=z2xgrid*hxgrid
  f2xobs=z2xobs*hxobs

  call InTrapCum(f2xgrid,xdelta,ngrid,f1xgrid)
  call IntFobs(f2xobs,f2xgrid,f1xgrid,xinxgrid,xidelta,nobs,ngrid,f1xobs)

  f1min=min(0.d0,minval(f1xgrid))

  call InTrapCum(f1xgrid,xdelta,ngrid,fxgrid)
  call IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid,fxobs)

  call IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid,fint)
  fxgrid=fxgrid-fint/xrange
  fxobs=fxobs-fint/xrange

  fxgrid=fpm*(fxgrid)+(alpha-f1min)*(xgrid-xmid)
  fxobs=fpm*(fxobs)+(alpha-f1min)*(xobs-xmid)

  f1xgrid=fpm*(f1xgrid)+alpha-f1min
  f1xobs=fpm*(f1xobs)+alpha-f1min

  f2xgrid=fpm*f2xgrid
  f2xobs=fpm*f2xobs

  return
end subroutine GetSf


subroutine GetRotateSf(fpm,omega,psi,alpha,theta,xobs,xgrid,phixobs,phixgrid,xdelta, &
                       xinxgrid,xidelta,xrange,xmid,quadfacts,intsimpfacts, &
                       nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: xinxgrid(nobs),quadfacts(nr,3),intsimpfacts(ngrid)
  real(8), intent(in) :: fpm,omega,psi,alpha,theta(nbasis+1),xobs(nobs),xgrid(ngrid)
  real(8), intent(in) :: phixobs(nr,nobs),phixgrid(nr,ngrid),xdelta,xidelta(nobs)
  real(8), intent(in) :: xrange,xmid

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  !internal arguments
  real(8) :: z2xobs(nobs),z2xgrid(ngrid),hxgrid(ngrid),hxobs(nobs)
  real(8) :: f2xgrid(ngrid),f2xobs(nobs),f1xgrid(ngrid),f1xobs(nobs),f1min,fint

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,z2xobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,z2xgrid)

  call SquishUp(xgrid,psi,omega,ngrid,hxgrid)
  call SquishUp(xobs,psi,omega,nobs,hxobs)

  f2xgrid=z2xgrid*hxgrid
  f2xobs=z2xobs*hxobs

  call InTrapCum(f2xgrid,xdelta,ngrid,f1xgrid)
  call IntFobs(f2xobs,f2xgrid,f1xgrid,xinxgrid,xidelta,nobs,ngrid,f1xobs)

  f1min=min(0.d0,minval(f1xgrid))

  call InTrapCum(f1xgrid,xdelta,ngrid,fxgrid)
  call IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid,fxobs)

  call IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid,fint)
  fxgrid=fxgrid-fint/xrange
  fxobs=fxobs-fint/xrange

  fxgrid=fpm*(fxgrid)+(alpha-f1min)*(xgrid-xmid)
  fxobs=fpm*(fxobs)+(alpha-f1min)*(xobs-xmid)

  f1xgrid=fpm*(f1xgrid)+alpha-f1min
  f1xobs=fpm*(f1xobs)+alpha-f1min

  f2xgrid=fpm*f2xgrid
  f2xobs=fpm*f2xobs

  return
end subroutine GetRotateSf


subroutine GetUf(fpm,omega,psi,theta,xobs,xgrid,phixobs,phixgrid,xdelta,xinxgrid, &
                 xidelta,xrange,quadfacts,intsimpfacts,nbasis,nr,nobs,ngrid,fxobs,fxgrid)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,nr,nobs,ngrid
  integer,intent(in) :: xinxgrid(nobs),quadfacts(nr,3),intsimpfacts(ngrid)
  real(8), intent(in) :: fpm,omega,psi,theta(nbasis+1),xobs(nobs),xgrid(ngrid)
  real(8), intent(in) :: phixobs(nr,nobs),phixgrid(nr,ngrid),xdelta,xidelta(nobs)
  real(8), intent(in) :: xrange

  !output arguments
  real(8),intent(out) :: fxobs(nobs),fxgrid(ngrid)

  !internal arguments
  real(8) :: z2xobs(nobs),z2xgrid(ngrid),hxgrid(ngrid),hxobs(nobs)
  real(8) :: f1xgrid(ngrid),f1xobs(nobs),fint

  call QuadMult(theta,phixobs,quadfacts,nbasis+1,nr,nobs,z2xobs)
  call QuadMult(theta,phixgrid,quadfacts,nbasis+1,nr,ngrid,z2xgrid)

  call SquishDown(xgrid,psi,omega,ngrid,hxgrid)
  call SquishDown(xobs,psi,omega,nobs,hxobs)

  f1xgrid=z2xgrid*hxgrid
  f1xobs=z2xobs*hxobs

  call InTrapCum(f1xgrid,xdelta,ngrid,fxgrid)
  call IntFobs(f1xobs,f1xgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid,fxobs)

  call IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid,fint)
  fxgrid=fxgrid-fint/xrange
  fxobs=fxobs-fint/xrange

  if (fpm.lt.0.d0) then
    fxgrid=-(fxgrid)
    fxobs=-(fxobs)

    f1xgrid=-(f1xgrid)
    f1xobs=-(f1xobs)

  end if

  return
end subroutine GetUf


subroutine SquishDown(x,psi,omega,n,xout)
  implicit none

  !input arguments
  integer,intent(in) :: n
  real(8), intent(in) :: x(n),psi,omega

  !output argument
  real(8),intent(out) :: xout(n)

  !internal arguments
  integer :: i
  real(8) :: c(n)

  c = psi*(x-omega)
  do i=1,n
    if(c(i).le.-100.d0) c(i)=-100.d0
    if(c(i).ge.100.d0) c(i)=100.d0
  end do
  c=dexp(c)
  xout=(1.d0-c)/(1.d0+c)

  return
end subroutine SquishDown


subroutine SquishUp(x,psi,omega,n,xout)
  implicit none

  !input arguments
  integer,intent(in) :: n
  real(8), intent(in) :: x(n),psi,omega

  !output argument
  real(8),intent(out) :: xout(n)

  !internal arguments
  integer :: i
  real(8) :: c(n)

  c = psi*(x-omega)
  do i=1,n
    if(c(i).le.-100.d0) c(i)=-100.d0
    if(c(i).ge.100.d0) c(i)=100.d0
  end do
  c=dexp(c)
  xout=(c-1.d0)/(c+1.d0)

  return
end subroutine SquishUp


subroutine InTrapCum(f,delta,n,fint)
  implicit none

  !input arguments
  integer,intent(in) :: n
  real(8), intent(in) :: f(n),delta

  !output argument
  real(8),intent(out) :: fint(n)

  !internal argument
  integer :: i
  real(8) :: fintmp(n)

  fintmp(1)=0.d0
  fintmp(2:n)=delta*(f(1:(n-1))+f(2:n))/2.d0

  fint=0.d0
  do i=1,n
    fint(i)=fint(i)+sum(fintmp(1:i))
  end do

  return
end subroutine InTrapCum



subroutine IntFobs(hobs,hxgrid,fxgrid,xinxgrid,xidelta,nobs,ngrid,fxobsout)
  implicit none

  !input arguments
  integer,intent(in) :: nobs,ngrid
  integer,intent(in) :: xinxgrid(nobs)
  real(8), intent(in) :: hobs(nobs),hxgrid(ngrid),fxgrid(ngrid),xidelta(nobs)

  !output argument
  real(8),intent(out) :: fxobsout(nobs)

  fxobsout=fxgrid(xinxgrid)
  fxobsout=fxobsout+xidelta*(hxgrid(xinxgrid)+hobs)/2.d0

  return
end subroutine IntFobs


subroutine IntSimpsonfxgrid(fxgrid,xdelta,intsimpfacts,ngrid,fint)
  implicit none

  !intput arguments
  integer,intent(in) :: ngrid
  integer,intent(in) :: intsimpfacts(ngrid)
  real(8), intent(in) :: fxgrid(ngrid),xdelta

  !output argument
  real(8),intent(out) :: fint

  fint=sum(fxgrid*dble(intsimpfacts))*xdelta/3.d0

  return
end subroutine IntSimpsonfxgrid


subroutine ConstFun(x,xrange,n,xout)
  implicit none

  !input arguments
  integer,intent(in) :: n
  real(8), intent(in) :: x(n),xrange

  !output argument
  real(8),intent(out) :: xout(n)

  xout=x*0.d0
  xout=1.d0/dsqrt(xrange)

  return
end subroutine ConstFun


subroutine ConstFun2(x,xmin,xrange,xout)
  implicit none

  !input arguments
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  xout=(x-xmin)*0.d0
  xout=1/xrange

  return
end subroutine ConstFun2


subroutine CosFun(x,xmin,xrange,n,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: n,nbasis
  real(8), intent(in) :: x(n),xmin,xrange

  !output argument
  real(8),intent(out) :: xout(n,nbasis)

  !internal arguments
  integer :: k
  real(8) :: z(n)

  z=(x-xmin)/xrange
  do k=1,nbasis
    xout(:,k)=dsqrt(2.d0/xrange)*dcos(dble(k)*PI*z)
  end do

  return
end subroutine CosFun


subroutine CosFun2(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=(2.d0/xrange)*(dcos(dble(kall)*PI*z)**2.d0)

  return
end subroutine CosFun2


subroutine ConstCosFun(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=(dsqrt(2.d0)/xrange)*dcos(dble(kall)*PI*z)

  return
end subroutine ConstCosFun


subroutine CrossProdFun(x,j,k,xmin,xrange,xout)
  implicit none

  !input arguments
  integer,intent(in) :: j,k
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=(2.d0/xrange)*dcos(dble(j)*PI*z)*dcos(dble(k)*PI*z)

  return
end subroutine CrossProdFun


subroutine IntConst2(x,xmin,xrange,xout)
  implicit none

  !input arguments
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  xout=(x-xmin)/xrange-1.d0/2.d0

  return
end subroutine IntConst2


subroutine IntIntConst2(x,xmin,xrange,xout)
  implicit none

  !input arguments
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  xout=(x-xmin)**2.d0/(2.d0*xrange)-xrange/6.d0

  return
end subroutine IntIntConst2


subroutine IntCos(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=dsqrt(2.d0)*dsin(dble(kall)*PI*z)/(PI*dble(kall))- &
       dsqrt(2.d0)*(1.d0-dcos(PI*dble(kall)))/((PI*dble(kall))**2.d0)

  return
end subroutine IntCos


subroutine IntIntCos(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal argument
  real(8) :: z

  z=(x-xmin)/xrange
  xout=dsqrt(2.d0)*xrange*(-dcos(dble(kall)*PI*z))/((PI*dble(kall))**2.d0)

  return
end subroutine IntIntCos


subroutine IntCos2(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=dsin(dble(kall)*2.d0*PI*z)/(2.d0*PI*dble(kall))+z-1.d0/2.d0

  return
end subroutine IntCos2


subroutine IntIntCos2(x,kall,xmin,xrange,nbasis,xout)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis
  integer,intent(in) :: kall(nbasis)
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout(nbasis)

  !internal argument
  real(8) :: z

  z=(x-xmin)/xrange
  xout=(1.d0-dcos(dble(kall)*2.d0*PI*z))*xrange/((2.d0*PI*dble(kall))**2.d0)+ &
       ((x-xmin)**2.d0)/(2.d0*xrange)-xrange/((2.d0*PI*dble(kall))**2.d0)-xrange/6.d0

  return
end subroutine IntIntCos2


subroutine IntCosCrossProd(x,j,k,xmin,xrange,xout)
  implicit none

  !input arguments
  integer,intent(in) :: j,k
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  !internal arguments
  real(8) :: z

  z=(x-xmin)/xrange
  xout=dsin(dble(j+k)*PI*z)/(PI*dble(j+k))+dsin(dble(k-j)*PI*z)/(PI*dble(k-j))- &
       (1.d0-dcos(PI*dble(j+k)))/((PI*dble(j+k))**2.d0)- &
       (1.d0-dcos(PI*dble(j-k)))/((PI*dble(j-k))**2.d0)

  return
end subroutine IntCosCrossProd


subroutine IntIntCrossProd(x,j,k,xmin,xrange,xout)
  implicit none

  !input arguments
  integer,intent(in) :: j,k
  real(8), intent(in) :: x,xmin,xrange

  !output argument
  real(8),intent(out) :: xout

  !internal argument
  real(8) :: z

  z=(x-xmin)/xrange
  xout=(1.d0-dcos(dble(j+k)*PI*z))*xrange/((PI*dble(j+k))**2.d0)+ &
       (1.d0-dcos(dble(j-k)*PI*z))*xrange/((PI*dble(j-k))**2.d0)- &
       xrange/((PI*dble(j+k))**2.d0)-xrange/((PI*dble(j-k))**2.d0)

  return
end subroutine IntIntCrossProd



subroutine GetPhi(x,xmin,xrange,phijj,phijk,phi00,phi0k,nbasis,n,phix)
  implicit none

  !input arguments
  integer,intent(in) :: nbasis,n
  real(8), intent(in) :: x(n),xmin,xrange
  external :: phijj,phijk,phi00,phi0k

  !output argument
  real(8),intent(out) :: phix((nbasis+1)*(nbasis+2)/2,n)

  !internal arguments
  integer :: kall(nbasis),nr,i,j,k
  real(8) :: xi,phi(nbasis,nbasis),phi1(nbasis+1,nbasis+1),a(nbasis),b,c,d(nbasis)

  phix=0.d0
  kall = (/ (k,k=1,nbasis) /)
  nr=(nbasis+1)*(nbasis+2)/2
  do i=1,n
    xi=x(i)
    call phijj(xi,kall,xmin,xrange,nbasis,a)
    call diagvec(a,nbasis,phi)
    do j=1,(nbasis-1)
      do k=(j+1),nbasis
        call phijk(xi,k,j,xmin,xrange,b)
        phi(j,k)=b
        phi(k,j)=b
      end do
    end do
    call phi00(xi,xmin,xrange,c)
    call phi0k(xi,kall,xmin,xrange,nbasis,d)
    phi1(1,1)=c
    phi1(1,2:nbasis)=d
    phi1(2:nbasis,1)=d
    phi1(2:nbasis,2:nbasis)=phi
    call vech(phi1,nbasis+1,nbasis+1,phix(:,i))
  end do

  return
end subroutine GetPhi

end module bsamTools
