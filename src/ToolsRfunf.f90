module ToolsRfunf
implicit none
real(8),parameter :: PI=3.141592653589793238462643383279502884197d0
contains

!*****************************************************************************************
! Random number gernerators
!*****************************************************************************************

function discrnd(n, probs)
  implicit none

  ! Input arguments
  integer, intent(in) :: n
  real(8), intent(in) :: probs(n)

  ! Output argument
  integer :: discrnd

  ! Internal argument
  real(8) :: rndunif,cum_probs(n),u
  integer :: i

  cum_probs=0.d0
  cum_probs(1)=probs(1)
  do i=2,n
    cum_probs(i)=cum_probs(i-1)+probs(i)
  end do

  u=rndunif()

  discrnd=1
  do i=1,n-1
    if (u .gt. cum_probs(i)) then
      discrnd = discrnd+1
    else
      discrnd = discrnd
      exit
    end if
  end do

  return
end function discrnd


subroutine mvnrnd(mu,cov,p,rn)
  implicit none

  ! Input arguments:
  integer, intent(in) :: p
  real(8), intent(in)  :: mu(p),cov(p,p)

  ! Output arguments:
  real(8), intent(out) :: rn(p)

  ! R functions:
  real(8) :: rndnorm

  ! Internal arguments:
  real(8) :: L(p,p),z(p)
  integer :: i,j,ok

  L=cov
  call dpotrf('L',p,L,p,ok)

  do i=1,p
    z(i)=rndnorm()
  enddo

  do i=1,p
    rn(i)=mu(i)
    do j=1,i
      rn(i)=rn(i)+L(i,j)*z(j)
    enddo
  enddo

  return
end subroutine mvnrnd


function rightmost_interval(u,lambda)
implicit none

!input arguments
real(8),intent(in) :: u,lambda

!output
integer :: rightmost_interval

!internal arguments
integer :: ok
real(8) :: j,R,Q

ok=0
R=1.d0
Q=dexp(-0.5d0*lambda)
j=0.d0
do
  call rchkusr() ! check interrupt

  j=j+1.d0
  R=R-((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
  if(R.gt.u) then
    ok=1
    exit
  end if
  j=j+1.d0
  R=R+((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
  if(R.lt.u) then
    ok=0
    exit
  end if
end do
rightmost_interval=ok

return
end function rightmost_interval


function leftmost_interval(u,lambda)
implicit none

!input arguments
real(8),intent(in) :: u,lambda

!output
integer :: leftmost_interval

!internal arguments
integer :: ok
real(8) :: j,R,Q,H,lu,K

ok=0
H=0.5d0*dlog(2.d0)+2.5d0*dlog(PI)-2.5d0*dlog(lambda)-&
  (PI**2.d0)/(2.d0*lambda)+0.5d0*lambda
lu=dlog(u)
R=1.d0
Q=dexp(-(PI**2.d0)/(2.d0*lambda))
K=lambda/(PI**2.d0)
j=0.d0
do
  call rchkusr() ! check interrupt

  j=j+1.d0
  R=R-K*(Q**((j**2.d0)-1.d0))
  if(H+dlog(R).gt.lu) then
    ok=1
    exit
  end if
  j=j+1.d0
  R=R+((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
  if(H+dlog(R).lt.lu) then
    ok=0
    exit
  end if
end do
leftmost_interval=ok

return
end function leftmost_interval



function gigrnd(lambda,psi,chi)
implicit none

! Input arguments
real(8),intent(in) :: lambda,psi,chi

! Output argument
real(8) :: gigrnd

! Internal arguments
real(8) :: gamrnd,invgaussrnd
real(8) :: a,b

! check parameters
if (chi.lt.0.d0) call rexit('chi must be non-negative')
if (psi.lt.0.d0) call rexit('psi must be non-negative')
if (chi.eq.0.d0) then
  if (lambda.le.0.d0) call rexit('lambda must be positive when chi = 0')
  if (psi.eq.0.d0) call rexit('psi and chi cannot both be 0')
end if
if (psi.eq.0.d0) then
  if (lambda.ge.0.d0) call rexit('lambda must be negative when psi = 0')
  if (chi.eq.0.d0) call rexit('chi and psi cannot both be 0')
end if

if (chi.eq.0.d0) then
  a=lambda
  b=psi/2.d0
  gigrnd=gamrnd(a,1.d0/b)
else if (psi.eq.0.d0) then
  a=-lambda
  b=chi/2.d0
  gigrnd=1.d0/gamrnd(a,1.d0/b)
else if (lambda.eq.-0.5d0) then
  b=chi
  a=dsqrt(b/psi)
  gigrnd=invgaussrnd(a,b)
else if (lambda.eq.1.d0) then
  gigrnd=rgig1(psi,chi)
else
  gigrnd=rgig(lambda,psi,chi)
end if

return
end function gigrnd


function rgig(lambda,psi,chi)
implicit none

! Input arguments
real(8),intent(in) :: lambda,psi,chi

! Output argument
real(8) :: rgig

! Internal arguments
real(8) :: tol,param(3),alpha,beta,m,upper,yM,yP,a,b,c,output,R1,R2,Y,Yc
real(8) :: rndunif,powerxy

tol=epsilon(lambda)

alpha=dsqrt(psi/chi)
beta=dsqrt(psi*chi)
m=(lambda-1.d0+dsqrt(powerxy(lambda-1.d0,2.d0)+powerxy(beta,2.d0)))/beta

param(1)=lambda
param(2)=beta
param(3)=m

upper=m
do
  if (gf(upper,param).gt.0.d0) exit
  upper=2*upper
end do

yM=zeroin(0.d0,m,param,tol)
yP=zeroin(m,upper,param,tol)

a=(yP-m)*powerxy(yP/m,0.5d0*(lambda-1.d0))*dexp(-0.25d0*beta*(yP+1.d0/yP-m-1.d0/m))
b=(yM-m)*powerxy(YM/m,0.5d0*(lambda-1.d0))*dexp(-0.25d0*beta*(yM+1.d0/yM-m-1.d0/m))
c=-0.25d0*beta*(m+1.d0/m)+0.5d0*(lambda-1.d0)*dlog(m)

output=0.d0
do
  R1=rndunif()
  R2=rndunif()
  Y=m+a*R2/R1+b*(1.d0-R2)/R1
  Yc=-0.5d0*(lambda-1.d0)*dlog(Y)+0.25d0*beta*(Y+1.d0/Y)+c
  if (Y.gt.0.d0 .and. -dlog(R1).ge.Yc) then
    output=Y
    exit
  end if
end do

rgig=output/alpha

return
end function rgig


function rgig1(psi,chi)
implicit none

! Input arguments
real(8),intent(in) :: psi,chi

! Output argument
real(8) :: rgig1

! Internal arguments
real(8) :: tol,param(3),alpha,beta,m,upper,yM,yP,a,b,c,output,R1,R2,Y
real(8) :: rndunif

tol=epsilon(psi)

alpha=dsqrt(psi/chi)
beta=dsqrt(psi*chi)
m=dabs(beta)/beta

param(1)=1.d0
param(2)=beta
param(3)=m

upper=m
do
  if (gf(upper,param).gt.0.d0) exit
  upper=2*upper
end do

yM=zeroin(0.d0,m,param,tol)
yP=zeroin(m,upper,param,tol)

a=(yP-m)*dexp(-0.25d0*beta*(yP+1.d0/yP-m-1.d0/m))
b=(yM-m)*dexp(-0.25d0*beta*(yM+1.d0/yM-m-1.d0/m))
c=-0.25d0*beta*(m+1.d0/m)

output=0.d0
do
  R1=rndunif()
  R2=rndunif()
  Y=m+a*R2/R1+b*(1.d0-R2)/R1
  if (Y.gt.0.d0 .and. -dlog(R1).ge.(0.25d0*beta*(Y+1.d0/Y)+c)) then
    output=Y
    exit
  end if
end do

rgig1=output/alpha

return
end function rgig1


function gf(x,param)
implicit none

! Input arguments
real(8),intent(in) :: x,param(3)

! Output argument
real(8) :: gf

! Internal arguments
real(8) :: powerxy
real(8) :: lambda,beta,m

lambda=param(1)
beta=param(2)
m=param(3)

gf=0.5d0*beta*powerxy(x,3.d0)-powerxy(x,2.d0)*(0.5d0*beta*m+lambda+1)+ &
   x*((lambda-1.d0)*m-0.5d0*beta)+0.5d0*beta*m

return
end function gf


!*****************************************************************************************
! Probability density function (pdf)
!*****************************************************************************************

function mvnpdf(x,mu,cov,d,log_p)
implicit none

! Input arguments
logical,intent(in) :: log_p
integer,intent(in) :: d
real(8),intent(in) :: x(d),mu(d),cov(d,d)

! Output argument
real(8) :: mvnpdf

! Internal arguments
real(8) :: detcov_half,covi(d,d),resid(d),logpdf,logconst
integer :: i,j,ok

covi=cov
call dpotrf('U',d,covi,d,ok)

detcov_half=1.d0
do i=1,d
  detcov_half=detcov_half*covi(i,i)
end do

call dpotri('U',d,covi,d,ok)
do i=1,(d-1)
  do j=(i+1),d
    covi(j,i)=covi(i,j)
  end do
end do

logconst=-dble(d)*dlog(2.d0*PI)/2.d0-dlog(detcov_half)
resid=x-mu
logpdf=logconst-dot_product(resid,matmul(covi,resid))/2.d0

mvnpdf=dexp(logpdf)
if (log_p) mvnpdf=logpdf

return
end function mvnpdf


!*****************************************************************************************
! Linear algebra
!*****************************************************************************************

subroutine inverse(R,p,Ri)
  implicit none

  ! Input arguments:
  integer, intent(in) :: p
  real(8), intent(in)  :: R(p,p)

  ! Output argument
  real(8), intent(out) :: Ri(p,p)

  ! Internal arguments:
  integer :: i,j,ok

  Ri=R
  call dpotrf('U',p,Ri,p,ok)
  call dpotri('U',p,Ri,p,ok)
  do i=1,(p-1)
    do j=(i+1),p
      Ri(j,i) = Ri(i,j)
    end do
  end do

  return
end subroutine inverse


function determinant(R,p)
  implicit none

  ! Input arguments
  integer, intent(in) :: p
  real(8), intent(in)  :: R(p,p)

  ! Output argument
  real(8) :: determinant

  ! Internal arguments
  integer :: i,ipiv(p),info
  real(8) :: Rlu(p,p)

  Rlu=R
  call dgetrf(p,p,Rlu,p,ipiv,info)

  determinant = 0.d0
  if (info.ne.0) then
    return
  endif
  determinant = 1.d0
  do i=1,p
    if (ipiv(i).ne.i) then
      determinant=-determinant*Rlu(i,i)
    else
      determinant=determinant*Rlu(i,i)
    endif
  end do

  return
end function determinant


subroutine covariance(A,n,p,cov)
implicit none

!input arguments
integer,intent(in) :: n,p
real(8), intent(in) :: A(n,p)

!output argument
real(8),intent(out) :: cov(p,p)

!internal arguments
real(8) :: z(n,p)

call sweep(A,n,p,.true.,.false.,z)
cov=matmul(transpose(z),z)/dble(n-1)

return
end subroutine covariance



subroutine sweep(x,n,p,center,scale,z)
implicit none

! input arguments
logical,intent(in) :: center,scale
integer,intent(in) :: n,p
real(8), intent(in) :: x(n,p)

! output argument
real(8),intent(out) :: z(n,p)

! internal arguments
real(8) :: xmean(p),xsd(p)
integer :: i

z=x
xsd=1.d0
xmean=0.d0

if (center) xmean=sum(x,1)/dble(n)
if (scale) xsd=(/ (dsqrt(sum((x(:,i)-xmean(i))**2.d0,1)/(dble(n)-1.d0)), i=1,n) /)

do i=1,p
  z(:,i)=(x(:,i)-xmean(i))/xsd(i)
end do

return
end subroutine sweep


subroutine Ikron(A,nra,nca,B,nrb,ncb,K)
  implicit none

  !input arguments
  integer,intent(in) :: nra,nca,nrb,ncb
  integer,intent(in) :: A(nra,nca),B(nrb,ncb)

  !output argument
  integer,intent(out) :: K(nra*nrb,nca*ncb)

  !internal arguments
  integer :: i,j

  K=0

  forall (i = 1:nra , j = 1:nca)
    K((nrb*(i-1)+1):(nrb*i),(ncb*(j-1)+1):(ncb*j))=A(i,j)*B
  end forall

  return
end subroutine Ikron


subroutine vech(mat,nr,nc,vec)
  implicit none

  !input arguments
  integer,intent(in) :: nr,nc
  real(8), intent(in) :: mat(nr,nc)

  !output argument
  real(8),intent(out) :: vec(nr * (nc + 1)/2)

  !internal arguemts
  integer :: ir,ic,k

  vec=0.d0

  k=1
  do ir=1,nr
    do ic=1,ir
      vec(k)=mat(ir,ic)
      k=k+1
    end do
  end do

  return
end subroutine vech


subroutine Ivech(mat,nr,nc,vec)
  implicit none

  !input arguments
  integer,intent(in) :: nr,nc
  integer,intent(in) :: mat(nr,nc)

  !output argument
  integer,intent(out) :: vec(nr * (nc + 1)/2)

  !internal arguemts
  integer :: ir,ic,k

  vec=0

  k=1
  do ir=1,nr
    do ic=1,ir
      vec(k)=mat(ir,ic)
      k=k+1
    end do
  end do

  return
end subroutine Ivech


subroutine diag(x,n,A)
  implicit none

  !input arguments
  real(8), intent(in) :: x
  integer,intent(in) :: n

  !output arguments
  real(8),intent(out) :: A(n,n)

  !internal argument
  integer :: i

  A=0.d0
  do i=1,n
    A(i,i)=x
  end do

  return
end subroutine diag


subroutine Idiag(x,n,A)
  implicit none

  !input arguments
  integer,intent(in) :: x
  integer,intent(in) :: n

  !output arguments
  integer,intent(out) :: A(n,n)

  !internal argument
  integer :: i

  A=0
  do i=1,n
    A(i,i)=x
  end do

  return
end subroutine Idiag


subroutine diagvec(x,n,A)
  implicit none

  !input arguments
  integer,intent(in) :: n
  real(8),intent(in) :: x(n)

  !output arguments
  real(8),intent(out) :: A(n,n)

  !internal argument
  integer :: i

  A=0.d0
  do i=1,n
    A(i,i)=x(i)
  end do

  return
end subroutine diagvec


subroutine find_uniquei(x, n, uniq_x, k)
  implicit none

  ! Input arguments
  integer, intent(in) :: n
  integer, intent(in) :: x(n)

  ! Output arguments
  integer, intent(out) :: uniq_x(n),k

  ! Internal arguments
  integer :: i, j, si, m

  uniq_x=0
  m=1
  do
    if (x(m) /= 0) exit
    m=m+1
  end do
  uniq_x(1)=x(m)
  si=m+1
  k=1
  outer: do i=si,n
    do j=1,k
      if (uniq_x(j) .eq. x(i)) then
        cycle outer
      endif
    end do
    if (x(i) /= 0) then
      k=k+1
      uniq_x(k)=x(i)
    end if
  end do outer

  return
end subroutine find_uniquei


subroutine which(logic,n,ind,k)
  implicit none

  ! Input arguments
  integer,intent(in) :: n
  logical,intent(in) :: logic(n)

  ! Output arguments
  integer,intent(out) :: ind(n),k

  ! Internal arguments
  integer :: i,l,x(n),y(n)

  ind=0
  k=count(logic)
  x=transfer(logic, (/ (1, i = 1, n) /))
  y=x*(/ (i, i = 1, n) /)
  l=1
  do i=1,n
    if (y(i) .ne. 0) then
      ind(l)=y(i)
      l=l+1
    end if
  end do

  return
end subroutine which


subroutine intsim(f,delta,n,fint)
implicit none

! Input arguments
integer,intent(in) :: n
real(8), intent(in) :: f(n),delta

! Output argument
real(8),intent(out) :: fint

! Internal arguments
integer :: i
real(8) :: t(n)

if(n.eq.2*floor(dble(n)/2.d0)) then
  call rexit('ERROR: Even number of rows for Simpson integration')
else if(n.eq.3) then
  t(1)=1.d0
  t(2)=4.d0
  t(3)=1.d0
  fint=sum(t*f)*delta/3.d0
else
  t(1)=1.d0
  do i=2,n-3,2
    t(i)=4.d0
    t(i+1)=2.d0
  end do
  t(n-1)=4.d0
  t(n)=1.d0
  fint=sum(t*f)*delta/3.d0
end if

return
end subroutine intsim



function zeroin(ax,bx,param,tol)
implicit none

! Input arguments
real(8),intent(in) :: ax,bx,tol
real(8),intent(inout) :: param(*)

! Output argument
real(8) :: zeroin

! Internal arguments
real(8) :: eps,tol1,a,b,c,d,e,fa,fb,fc,xm,p,q,r,s

eps = 1.d0
10 eps = eps/2.d0
tol1 = 1.d0 + eps
if (tol1 .gt. 1.d0) go to 10

a = ax
b = bx
fa = gf(a,param)
fb = gf(b,param)

20 c = a
fc = fa
d = b - a
e = d
30 if (dabs(fc) .ge. dabs(fb)) go to 40
a = b
b = c
c = a
fa = fb
fb = fc
fc = fa

40 tol1 = 2.d0*eps*dabs(b) + 0.5*tol
xm = .5*(c - b)
if (dabs(xm) .le. tol1) go to 90
if (fb .eq. 0.d0) go to 90

if (dabs(e) .lt. tol1) go to 70
if (dabs(fa) .le. dabs(fb)) go to 70

if (a .ne. c) go to 50

s = fb/fa
p = 2.d0*xm*s
q = 1.d0 - s
go to 60

50 q = fa/fc
r = fb/fc
s = fb/fa
p = s*(2.d0*xm*q*(q - r) - (b - a)*(r - 1.d0))
q = (q - 1.d0)*(r - 1.d0)*(s - 1.d0)

60 if (p .gt. 0.d0) q = -q
p = dabs(p)

if ((2.d0*p) .ge. (3.d0*xm*q - dabs(tol1*q))) go to 70
if (p .ge. dabs(0.5*e*q)) go to 70
e = d
d = p/q
go to 80

70 d = xm
e = d

80 a = b
fa = fb
if (dabs(d) .gt. tol1) b = b + d
if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
fb = gf(b,param)
if ((fb*(fc/dabs(fc))) .gt. 0.d0) go to 20
go to 30

90 zeroin = b

return
end function zeroin

end module ToolsRfunf
