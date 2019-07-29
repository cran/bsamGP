subroutine bsad(verbose,nblow,smcmc,nskip,ndisp,kappaloop,nobs,nint,npar,MaxNCos,&
                probk,cdfk,np,gmax,smooth,PriorProbs,nmodels,&
                BetaVari0,BetaVariMean0,r0,s0,u0,v0,dmat,ddata,dtd,&
                phi,phidata,ndata,xdelta,lnPriorKappa,cdfPriorKappa,KappaGrid,&
                betaParg,sigmaParg,betag,sigmag,taug,gamg,thetag,kappag,&
                betaMaxKappag,sigmaMaxKappag,tauMaxKappag,gamMaxKappag,&
                thetaMaxKappag,fparg,fsemig,fsemiMaxKappag,Probzzg,ExactLogLikeg,&
                MetProbPar,MetProbSemi)
use ToolsRfunf
implicit none

integer,intent(in) :: verbose,nblow,smcmc,nskip,ndisp,nobs,npar,nint,MaxNCos,smooth
integer,intent(in) :: kappaloop,np,nmodels,ndata(nint),KappaGrid(MaxNCos+1)
real(8), intent(in) :: probk(np,np),cdfk(np,np),gmax
real(8), intent(in) :: PriorProbs(nmodels),BetaVari0(npar,npar)
real(8), intent(in) :: BetaVariMean0(npar),r0,s0,u0,v0
real(8), intent(in) :: dmat(nint,npar),ddata(nobs,npar),xdelta,phi(nint,MaxNCos)
real(8), intent(in) :: phidata(nobs,MaxNCos),dtd(npar,npar),lnPriorKappa(MaxNCos+1)
real(8), intent(in) :: cdfPriorKappa(MaxNCos+1)

integer,intent(out) :: MetProbPar,MetProbSemi,kappag(smcmc)
real(8), intent(out) :: betaParg(smcmc,npar),sigmaParg(smcmc),betag(smcmc,npar)
real(8), intent(out) :: taug(smcmc),gamg(smcmc),thetag(smcmc,MaxNCos)
real(8), intent(out) :: betaMaxKappag(smcmc,npar),sigmaMaxKappag(smcmc),sigmag(smcmc)
real(8), intent(out) :: tauMaxKappag(smcmc),gamMaxKappag(smcmc),ExactLogLikeg(smcmc,3)
real(8), intent(out) :: thetaMaxKappag(smcmc,MaxNCos),fparg(smcmc,nint)
real(8), intent(out) :: fsemig(smcmc,nint),fsemiMaxKappag(smcmc,nint)
real(8), intent(out) :: Probzzg(smcmc,nmodels)

real(8)  :: stime,itime
integer :: ikappa,imcmc,isave,nmcmc

integer :: k
real(8) :: ck,w0,wk,kall(MaxNCos)

integer :: kappa
real(8) :: betaPar(npar),sigmaPar,sigma2Par,ydataPar(nint),loglike,dmatt(npar,nint)
real(8) :: beta(npar),sigma,sigma2,tau2,gam,theta(MaxNCos),ydata(nint)
real(8) :: betaMaxKappa(npar),sigmaMaxKappa,sigma2MaxKappa,tau2MaxKappa
real(8) :: gamMaxKappa,thetaMaxKappa(MaxNCos),ydataMaxKappa(nint),fx(nint),fint

real(8) :: ymeanPar(nint),c,ExactLogLikeParNew,ExactLogLikeParOld,yint
real(8) :: ymean(nint),ExactLogLikeSemi,ExactLogLikeSemiOld,ExactLogLikeSemiNew

real(8) :: ymeanMaxKappa(nint),ExactLogLikeSemiMaxKappa
real(8) :: ExactLogLikeSemiOldMaxKappa,ExactLogLikeSemiNewMaxKappa,Probzz(nmodels)

call cpu_time(stime)
call rndstart()

dmatt=transpose(dmat)

if(smooth.eq.1) then
  kall=(/ (dble(k),k=1,MaxNCos) /)
else
  kall=(/ (dlog(dble(k)),k=1,MaxNCos) /)
end if
ck=sum(kall)/2.d0 ;
w0=1.d0
wk=ck-w0

Probzz=0.5d0

betaPar=0.d0
sigmaPar=1.d0
sigma2Par=sigmaPar**2.d0
ydataPar=0.d0

beta=0.d0
sigma=1.d0
sigma2=sigma**2.d0
tau2=1.d0
gam=1.d0
theta=0.d0
ydata=0.d0
kappa=MaxNCos/2

betaMaxKappa=0.d0
sigmaMaxKappa=1.d0
sigma2MaxKappa=1.d0
tau2MaxKappa=1.d0
gamMaxKappa=1.d0
thetaMaxKappa=0.d0
ydataMaxKappa=0.d0

ymeanPar=matmul(dmat,betaPar)
c=maxval(ymeanPar)-3.d0
call intsim(dexp(ymeanPar-c),xdelta,nint,yint)
ExactLogLikeParNew=sum(matmul(ddata,betaPar))-dble(nobs)*(c+dlog(yint))
ExactLogLikeParOld=ExactLogLikeParNew

ymean=matmul(dmat,beta)
ExactLogLikeSemi=sum(matmul(ddata,beta))
if (kappa.gt.0) then
  ymean=ymean+matmul(phi(:,1:kappa),theta(1:kappa))
  ExactLogLikeSemi=ExactLogLikeSemi+sum(matmul(phidata(:,1:kappa),theta(1:kappa)))
end if
c=maxval(ymean)-3.d0
call intsim(dexp(ymean-c),xdelta,nint,yint)
ExactLogLikeSemi=ExactLogLikeSemi-dble(nobs)*(c+dlog(yint))
ExactLogLikeSemiOld=ExactLogLikeSemi
ExactLogLikeSemiNew=ExactLogLikeSemi

ymeanMaxKappa=matmul(dmat,betaMaxKappa)+matmul(phi,thetaMaxKappa)
ExactLogLikeSemiMaxKappa=sum(matmul(ddata,betaMaxKappa)+matmul(phidata,thetaMaxKappa))
c=maxval(ymeanMaxKappa)-3.d0
call intsim(dexp(ymeanMaxKappa-c),xdelta,nint,yint)
ExactLogLikeSemiMaxKappa=ExactLogLikeSemiMaxKappa-dble(nobs)*(c+dlog(yint))
ExactLogLikeSemiOldMaxKappa=ExactLogLikeSemiMaxKappa
ExactLogLikeSemiNewMaxKappa=ExactLogLikeSemiMaxKappa

if (verbose.eq.1) then
  call dblepr('Burnin ...',-1,1.d0,0)
end if
MetProbPar=0
MetProbSemi=0
do imcmc=1,nblow
  call rchkusr()

  call DoMCMCPar()
  call DoMCMCSemi()
  call DoMCMCSemiMaxKappa()
  do ikappa=1,kappaloop
    call DoMCMCPar()
    call DoMCMCSemiNoKappa()
    call DoMCMCSemiMaxKappa()
  end do
end do

if (verbose.eq.1) then
  call dblepr('Main iterations ...',-1,1.d0,0)
end if
MetProbPar=0
MetProbSemi=0
isave=1
nmcmc=nskip*smcmc
do imcmc=1,nmcmc
  call rchkusr()

  call DoMCMCPar()
  call DoMCMCSemi()
  call DoMCMCSemiMaxKappa()
  do ikappa=1,kappaloop
    call DoMCMCPar()
    call DoMCMCSemiNoKappa()
    call DoMCMCSemiMaxKappa()
  end do
  call ProbStay()

  if(mod(imcmc,nskip).eq.0) then
    ExactLogLikeg(isave,1)=ExactLogLikeParOld
    ExactLogLikeg(isave,2)=ExactLogLikeSemiOld
    ExactLogLikeg(isave,3)=ExactLogLikeSemiOldMaxKappa
    Probzzg(isave,:)=Probzz

    betaParg(isave,:)=betaPar
    sigmaParg(isave)=dsqrt(sigma2Par)
    ymeanPar=matmul(dmat,betaPar)
    fx=dexp(ymeanPar-maxval(ymeanPar))
    call intsim(fx,xdelta,nint,fint)
    fx=fx/fint
    fparg(isave,:)=fx

    betag(isave,:)=beta
    sigmag(isave)=dsqrt(sigma2)
    taug(isave)=dsqrt(tau2)
    gamg(isave)=gam
    kappag(isave)=kappa
    thetag(isave,:)=theta
    ymean=matmul(dmat,beta)
    if(kappa.gt.0) then
      ymean=ymean+matmul(phi(:,1:kappa),theta(1:kappa))
    end if
    fx=dexp(ymean-maxval(ymean))
    call intsim(fx,xdelta,nint,fint)
    fx=fx/fint
    fsemig(isave,:)=fx

    betaMaxKappag(isave,:)=betaMaxKappa
    sigmaMaxKappag(isave)=dsqrt(sigma2MaxKappa)
    tauMaxKappag(isave)=dsqrt(tau2MaxKappa)
    gamMaxKappag(isave)=gamMaxKappa
    thetaMaxKappag(isave,:)=thetaMaxKappa
    ymean=matmul(dmat,betaMaxKappa)+matmul(phi,thetaMaxKappa)
    fx=dexp(ymean-maxval(ymean))
    call intsim(fx,xdelta,nint,fint)
    fx=fx/fint
    fsemiMaxKappag(isave,:)=fx

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

contains

subroutine DoMCMCPar()
implicit none

real(8) :: betaParNew(npar),ymeanParNew(nint),c,yint
real(8) :: residOld(nint),residNew(nint),testp,rndunif

ymeanPar=matmul(dmat,betaPar)
call GetLogits(ndata,ymeanPar,sigma2Par,nobs,nint,ydataPar,loglike)

betaParNew=betaPar
call GetPar(ydataPar,dmat,dmatt,dtd,BetaVari0,BetaVariMean0,r0,s0,npar,nint, &
            betaParNew,sigma2Par)

ymeanParNew=matmul(dmat,betaParNew)
c=maxval(ymeanParNew)-3.d0
call intsim(dexp(ymeanParNew-c),xdelta,nint,yint)
ExactLogLikeParNew=sum(matmul(ddata,betaParNew))-dble(nobs)*(c+dlog(yint))

residOld=ydataPar-ymeanPar
residNew=ydataPar-ymeanParNew

testp=ExactLogLikeParNew-ExactLogLikeParOld- &
      sum(residOld**2.d0)/(2.d0*sigma2Par)+sum(residNew**2.d0)/(2.d0*sigma2Par)

if(dlog(rndunif()).lt.testp) then
  betaPar=betaParNew
  ExactLogLikeParOld=ExactLogLikeParNew
  ymeanPar=ymeanParNew
  MetProbPar=MetProbPar+1
end if

return
end subroutine DoMCMCPar


subroutine DoMCMCSemi()
implicit none

integer :: kappaNew
real(8) :: phitheta(nint),resid(nint),betaNew(npar),lnpkappa,thetaNew(MaxNCos)
real(8) :: ymeanNew(nint),phithetaNew(nint),c,yint,residOld(nint),residNew(nint)
real(8) :: residbotn(nint),residbnto(nint),testp,rndunif

ymean=matmul(dmat,beta)
phitheta=0.d0
if(kappa.gt.0) then
  phitheta=matmul(phi(:,1:kappa),theta(1:kappa))
  ymean=ymean+phitheta
end if
call GetLogits(ndata,ymean,sigma2,nobs,nint,ydata,loglike)

resid=ydata-phitheta
betaNew=beta
call GetPar(resid,dmat,dmatt,dtd,BetaVari0,BetaVariMean0,r0,s0,npar,nint,betaNew,sigma2)

call GetKappa(kappa,probk,cdfk,MaxNCos,np,kappaNew,lnpkappa)

resid=ydata-matmul(dmat,beta)
call GetFourier(resid,phi,kappaNew,sigma2,tau2,gam,kall,MaxNCos,nint,thetaNew)
ymeanNew=matmul(dmat,betaNew)

phithetaNew=0.d0
ExactLogLikeSemiNew=sum(matmul(ddata,betaNew))
if(kappaNew.gt.0) then
  phithetaNew=matmul(phi(:,1:kappaNew),thetaNew(1:kappaNew))
  ymeanNew=ymeanNew+phithetaNew
  ExactLogLikeSemiNew=ExactLogLikeSemiNew+sum(matmul(phidata(:,1:kappaNew),thetaNew(1:kappaNew)))
end if

c=maxval(ymeanNew)-3.d0
call intsim(dexp(ymeanNew-c),xdelta,nint,yint)
ExactLogLikeSemiNew=ExactLogLikeSemiNew-dble(nobs)*(c+dlog(yint))

residOld=ydata-ymean
residNew=ydata-ymeanNew
residbotn=ydata-matmul(dmat,beta)-phithetaNew
residbnto=ydata-matmul(dmat,betaNew)-phitheta

testp=ExactLogLikeSemiNew-ExactLogLikeSemiOld+ &
      lnPriorKappa(kappaNew+1)-lnPriorKappa(kappa+1)+lnpkappa-&
      sum(residOld**2.d0)/(2.d0*sigma2)-sum(residbotn**2.d0)/(2.d0*sigma2)+&
      sum(residNew**2.d0)/(2.d0*sigma2)+sum(residbnto**2.d0)/(2.d0*sigma2)

if(dlog(rndunif()).lt.testp) then
  beta=betaNew
  kappa=kappaNew
  theta=thetaNew
  ExactLogLikeSemiOld=ExactLogLikeSemiNew
  MetProbSemi=MetProbSemi+1
  ymean=ymeanNew
end if

call GetSmooth(theta,gmax,kall,wk,u0,v0,MaxNCos,tau2,gam)

return
end subroutine DoMCMCSemi


subroutine DoMCMCSemiNoKappa()
implicit none

real(8) :: phitheta(nint),resid(nint),betaNew(npar),thetaNew(MaxNCos)
real(8) :: ymeanNew(nint),phithetaNew(nint),c,yint,residOld(nint),residNew(nint)
real(8) :: residbotn(nint),residbnto(nint),testp,rndunif

ymean=matmul(dmat,beta)
phitheta=0.d0
if(kappa.gt.0) then
  phitheta=matmul(phi(:,1:kappa),theta(1:kappa))
  ymean=ymean+phitheta
end if
call GetLogits(ndata,ymean,sigma2,nobs,nint,ydata,loglike)

resid=ydata-phitheta
betaNew=beta
call GetPar(resid,dmat,dmatt,dtd,BetaVari0,BetaVariMean0,r0,s0,npar,nint,betaNew,sigma2)

resid=ydata-matmul(dmat,beta)
call GetFourier(resid,phi,kappa,sigma2,tau2,gam,kall,MaxNCos,nint,thetaNew)

ymeanNew=matmul(dmat,betaNew)
phithetaNew=0.d0
ExactLogLikeSemiNew=sum(matmul(ddata,betaNew))
if(kappa.gt.0) then
  phithetaNew=matmul(phi(:,1:kappa),thetaNew(1:kappa))
  ymeanNew=ymeanNew+phithetaNew
  ExactLogLikeSemiNew=ExactLogLikeSemiNew+sum(matmul(phidata(:,1:kappa),thetaNew(1:kappa)))
end if

c=maxval(ymeanNew)-3.d0
call intsim(dexp(ymeanNew-c),xdelta,nint,yint)
ExactLogLikeSemiNew=ExactLogLikeSemiNew-dble(nobs)*(c+dlog(yint))

residOld=ydata-ymean
residNew=ydata-ymeanNew
residbotn=ydata-matmul(dmat,beta)-phithetaNew
residbnto=ydata-matmul(dmat,betaNew)-phitheta

testp=ExactLogLikeSemiNew-ExactLogLikeSemiOld- &
      sum(residOld**2.d0)/(2.d0*sigma2)-sum(residbotn**2.d0)/(2.d0*sigma2)+&
      sum(residNew**2.d0)/(2.d0*sigma2)+sum(residbnto**2.d0)/(2.d0*sigma2)

if(dlog(rndunif()).lt.testp) then
  beta=betaNew
  theta=thetaNew
  ExactLogLikeSemiOld=ExactLogLikeSemiNew
  MetProbSemi=MetProbSemi+1
  ymean=ymeanNew
end if

call GetSmooth(theta,gmax,kall,wk,u0,v0,MaxNCos,tau2,gam)

return
end subroutine DoMCMCSemiNoKappa


subroutine DoMCMCSemiMaxKappa()
implicit none

real(8) :: phitheta(nint),resid(nint),betaNew(npar),thetaNew(MaxNCos)
real(8) :: ymeanNew(nint),phithetaNew(nint),c,yint,residOld(nint),residNew(nint)
real(8) :: residbotn(nint),residbnto(nint),testp,rndunif

phitheta=matmul(phi,thetaMaxKappa)
ymeanMaxKappa=matmul(dmat,betaMaxKappa)+phitheta
call GetLogits(ndata,ymeanMaxKappa,sigma2MaxKappa,nobs,nint,ydataMaxKappa,loglike)

resid=ydataMaxKappa-phitheta
betaNew=betaMaxKappa
call GetPar(resid,dmat,dmatt,dtd,BetaVari0,BetaVariMean0,r0,s0,npar,nint, &
            betaNew,sigma2MaxKappa)

resid=ydataMaxKappa-matmul(dmat,betaMaxKappa)
call GetFourier(resid,phi,MaxNCos,sigma2MaxKappa,tau2MaxKappa,gamMaxKappa,kall, &
                MaxNCos,nint,thetaNew)
phithetaNew=matmul(phi,thetaNew)
ymeanNew=matmul(dmat,betaNew)+phithetaNew

ExactLogLikeSemiNewMaxKappa=sum(matmul(ddata,betaNew))+sum(matmul(phidata,thetaNew))

c=maxval(ymeanNew)-3.d0
call intsim(dexp(ymeanNew-c),xdelta,nint,yint)
ExactLogLikeSemiNewMaxKappa=ExactLogLikeSemiNewMaxKappa-dble(nobs)*(c+dlog(yint))

residOld=ydataMaxKappa-ymeanMaxKappa
residNew=ydataMaxKappa-ymeanNew
residbotn=ydataMaxKappa-matmul(dmat,betaMaxKappa)-phithetaNew
residbnto=ydataMaxKappa-matmul(dmat,betaNew)-phitheta

testp=ExactLogLikeSemiNewMaxKappa-ExactLogLikeSemiOldMaxKappa- &
      sum(residOld**2.d0)/(2.d0*sigma2MaxKappa)- &
      sum(residbotn**2.d0)/(2.d0*sigma2MaxKappa)+ &
      sum(residNew**2.d0)/(2.d0*sigma2MaxKappa)+ &
      sum(residbnto**2.d0)/(2.d0*sigma2MaxKappa)

if(dlog(rndunif()).lt.testp) then
  betaMaxKappa=betaNew
  thetaMaxKappa=thetaNew
  ExactLogLikeSemiOldMaxKappa=ExactLogLikeSemiNewMaxKappa
  MetProbSemi=MetProbSemi+1
  ymeanMaxKappa=ymeanNew
end if

call GetSmooth(thetaMaxKappa,gmax,kall,wk,u0,v0,MaxNCos,tau2MaxKappa,gamMaxKappa)

return
end subroutine DoMCMCSemiMaxKappa


real(8) function GetLogLike(ndata,ydata,nint,nobs)
implicit none

integer,intent(in) :: nint,ndata(nint),nobs
real(8), intent(in) :: ydata(nint)

real(8) :: expy(nint)

expy=dexp(ydata)
GetLogLike=dot_product(dble(ndata),ydata)-dble(nobs)*dlog(sum(expy))

return
end function GetLogLike


subroutine GetLogits(ndata,ymean,sigma2,nobs,nint,ydata,loglike)
implicit none

integer,intent(in) :: nint,nobs
integer,intent(in) :: ndata(nint)
real(8), intent(in) :: ymean(nint),sigma2

real(8),intent(inout) :: ydata(nint)
real(8),intent(out) :: loglike

integer :: j,ymiss,ismiss
real(8) :: ymin,ymax,ymaxx,maxold,mydata,ym(nint),ys,expy(nint)
real(8) :: sexpy,uall(nint),u,v,expyold,yj,expynew,missy(3)
real(8) :: powerxy,rndunif,tnormrnd

ymin=-10.d0
ymaxx=10.d0
maxold=maxval(ydata)
mydata=sum(ydata)/dble(nint)
ym=ymean+sigma2*dble(ndata)
ys=dsqrt(sigma2)
expy=dexp(ydata)
sexpy=sum(expy)
uall=(/ (rndunif(),j=1,nint) /)
do j=1,nint
  u=uall(j)
  v=powerxy(u,(1.d0/dble(nobs)))
  expyold=expy(j)
  ymax=dlog(dabs(sexpy/v-(sexpy-expyold)))
  ymax=min(ymax,ymaxx)
  yj=tnormrnd(ym(j),dsqrt(sigma2),ymin,ymax)
  ymiss=ismiss(yj)
  if(ymiss.eq.1) then
    missy(1)=ydata(j)
    missy(2)=ymean(j)
    missy(3)=ymax
    call dblepr('Missing ydata: ',-1,missy,3)
    yj=mydata
  end if
  ydata(j)=yj
  expynew=dexp(yj)
  expy(j)=expynew
  sexpy=dabs(sexpy-expyold+expynew)
end do

loglike=dot_product(dble(ndata),ydata)-dble(nobs)*dlog(sexpy)

return
end subroutine GetLogits


subroutine GetSmooth(theta,gmax,kall,wk,u0,v0,MaxNCos,tau2,gam)
implicit none

integer,intent(in) :: MaxNCos
real(8), intent(in) :: theta(MaxNCos),gmax,kall(MaxNCos),wk,u0,v0

real(8),intent(inout) :: tau2,gam

integer :: i,j,k,nloop,NOTzero
real(8) :: un,gamvec(MaxNCos),ctol,gmin,theta2(MaxNCos)
real(8) :: th2gam(MaxNCos),vn,ueta,gmaxx,rndunif,gamrnd
real(8),allocatable :: th2gam1(:),kall1(:),ueta1(:),hylat(:),c2(:)


un=u0+dble(MaxNCos)
gamvec=dexp(-kall*gam)
ctol=0.00001d0

gmin=-dlog(ctol)/dble(MaxNCos)

nloop=50
theta2=dabs(theta)**2.d0

do j=1,nloop
  th2gam=theta2/gamvec
  vn=v0+sum(th2gam)
  tau2=vn/(2.d0*gamrnd(un/2.d0,1.d0))

  NOTzero=count(th2gam.ne.0.d0)
  if (NOTzero.eq.0) then
    ueta=rndunif()
    gam=gmax+dlog(ueta+(1.d0-ueta)*dexp(-wk*(gmax-gmin)))/wk
  else
    allocate(th2gam1(NOTzero),kall1(NOTzero),c2(NOTzero))
    allocate(ueta1(NOTzero),hylat(NOTzero+1))
    k=1
    do i=1,MaxNCos
      if(th2gam(i).ne.0.d0) then
        th2gam1(k)=th2gam(i)
        kall1(k)=kall(i)
        k=k+1
      end if
    end do
    ueta1=(/ (rndunif(),i=1,NOTzero) /)
    c2=(2.d0*tau2/th2gam1)
    hylat(1:NOTzero)=gam+dlog(1.d0-c2*dlog(ueta1))/kall1
    hylat(NOTzero+1)=gmax
    gmaxx=minval(hylat)

    ueta=rndunif()
    gam=gmaxx+dlog(ueta+(1.d0-ueta)*dexp(-wk*(gmaxx-gmin)))/wk
    deallocate(th2gam1,kall1,ueta1,hylat,c2)
  end if
end do

return
end subroutine GetSmooth


subroutine GetFourier(ydata,phi,kappa,sigma2,tau2,gam,kall,MaxNCos,nint,theta)
implicit none

integer,intent(in) :: kappa,MaxNCos,nint
real(8), intent(in) :: ydata(nint),phi(nint,MaxNCos),sigma2,tau2,gam,kall(MaxNCos)

real(8),intent(out) :: theta(MaxNCos)

integer :: i
real(8) :: gamvec(MaxNCos),zt(MaxNCos),phiy(MaxNCos)
real(8) :: v2i(MaxNCos),v2(MaxNCos),rndnorm

gamvec=dexp(-kall*gam)
if (kappa.eq.0) then
  theta=(/ (rndnorm(),i=1,MaxNCos) /)
  theta=dsqrt(tau2*gamvec)*theta
else
  zt=0.d0
  zt(1:kappa)=1.d0
  phiy(1:kappa)=matmul(transpose(phi(:,1:kappa)),ydata)/sigma2
  if (kappa.lt.MaxNCos) then
    phiy((kappa+1):MaxNCos)=0.d0
  end if
  v2i=(dble(nint)/sigma2)*zt+1.d0/(tau2*gamvec)
  v2=1.d0/v2i
  theta=(/ (rndnorm(),i=1,MaxNCos) /)
  theta=v2*phiy+dsqrt(v2)*theta
end if

return
end subroutine GetFourier


subroutine GetPar(ydata,dmat,dmatt,dtd,BetaVari0,BetaVariMean0,r0,s0,npar,nint,beta,sigma2)
implicit none

integer,intent(in) :: npar,nint
real(8), intent(in) :: ydata(nint),dmat(nint,npar),dmatt(npar,nint),dtd(npar,npar)
real(8), intent(in) :: BetaVari0(npar,npar),BetaVariMean0(npar),r0,s0

real(8),intent(inout) :: beta(npar)
real(8),intent(out) :: sigma2

real(8) :: rn,resid(nint),sn,gamrnd
real(8) :: Bni(npar,npar),Bn(npar,npar),bhat(npar)

rn=r0+dble(nint)
resid=dabs(ydata-matmul(dmat,beta))
sn=s0+sum(resid**2.d0)
sigma2=sn/(2.d0*gamrnd(rn/2.d0,1.d0))

Bni=dtd/sigma2+BetaVari0
if(determinant(Bni,npar).eq.0.d0) then
  sigma2=1.d0
  Bni=dtd/sigma2+BetaVari0
end if
call inverse(Bni,npar,Bn)
bhat=matmul(Bn,matmul(dmatt,ydata)/sigma2+BetaVariMean0)
call mvnrnd(bhat,Bn,npar,beta)

return
end subroutine GetPar


subroutine GetKappa(kappa,probk,cdfk,MaxNCos,np,KappaNew,lnpkappa)
implicit none

integer,intent(in) :: kappa,MaxNCos,np
real(8), intent(in) :: probk(np,np),cdfk(np,np)

integer,intent(out) :: KappaNew
real(8), intent(out) :: lnpkappa

integer :: i,z,bw,delta
real(8) :: u,pKappaNew,pKappa,rndunif
real(8), allocatable :: acdf(:),aprob(:)
integer,allocatable :: kgrid(:)

bw=floor(dble(np)/2.d0)
z=0
if(kappa.lt.bw) then
  allocate(acdf(kappa+1+bw),aprob(kappa+1+bw),kgrid(kappa+1+bw))
  acdf=cdfk(1:(kappa+1+bw),kappa+1)
  aprob=probk(1:(kappa+1+bw),kappa+1)
  kgrid=(/ (i,i=0,kappa+bw) /)
  u=rndunif()
  if(u.le.acdf(1)) then
    z=1
  else if(u.gt.acdf(kappa+1+bw)) then
    z=kappa+1+bw
  else
    do i=2,(kappa+1+bw)
      if(acdf(i-1).lt.u .and. u.le.acdf(i)) then
        z=i
        exit
      end if
    end do
  end if
  KappaNew=kgrid(z)
  pKappaNew=aprob(z)
  deallocate(acdf,aprob,kgrid)
else if(kappa.gt.(MaxNCos-bw)) then
  allocate(acdf(np-(1+bw+kappa-MaxNCos)+1),aprob(np-(1+bw+kappa-MaxNCos)+1))
  allocate(kgrid(np-(1+bw+kappa-MaxNCos)+1))
  acdf=cdfk((1+bw+kappa-MaxNCos):np,np-MaxNCos+kappa)
  aprob=probk((1+bw+kappa-MaxNCos):np,np-MaxNCos+kappa)
  kgrid=(/ (i,i=1,(np-(1+bw+kappa-MaxNCos)+1)) /)
  kgrid=kgrid+(kappa-bw-1)
  u=rndunif()
  if(u.le.acdf(1)) then
    z=1
  else if(u.gt.acdf(np-(1+bw+kappa-MaxNCos)+1)) then
    z=np-(1+bw+kappa-MaxNCos)+1
  else
    do i=2,np-(1+bw+kappa-MaxNCos)+1
      if(acdf(i-1).lt.u .and. u.le.acdf(i)) then
        z=i
        exit
      end if
    end do
  end if
  KappaNew=kgrid(z)
  pKappaNew=aprob(z)
  deallocate(acdf,aprob,kgrid)
else
  allocate(acdf(np),aprob(np),kgrid(np))
  acdf=cdfk(:,bw+1)
  aprob=probk(:,bw+1)
  kgrid=(/ (i,i=1,np) /)
  kgrid=kgrid+(kappa-bw-1)
  u=rndunif()
  if(u.le.acdf(1)) then
    z=1
  else if(u.gt.acdf(np)) then
    z=np
  else
    do i=2,np
      if(acdf(i-1).lt.u .and. u.le.acdf(i)) then
        z=i
        exit
      end if
    end do
  end if
  KappaNew=kgrid(z)
  pKappaNew=aprob(z)
  deallocate(acdf,aprob,kgrid)
end if

delta=iabs(KappaNew-kappa)
if(KappaNew.lt.bw) then
  pKappa=probk(KappaNew+1+delta,KappaNew+1)
else if(KappaNew.gt.(MaxNCos-bw)) then
  pKappa=probk(np-MaxNCos+KappaNew-delta,np-MaxNCos+KappaNew)
else
  pKappa=probk(bw+1+delta,bw+1)
end if
lnpkappa=dlog(pKappa)-dlog(pKappaNew)

return
end subroutine GetKappa


subroutine ProbStay()
implicit none

integer :: i,z,kappaPar
real(8) :: u,ctol,gamPar,tau2Par,tauPar,ymean(nint),loglike12,c,fint
real(8) :: loglike21,gmin
real(8) :: infp(nmodels),probz(nmodels)
real(8) :: rndunif,gamrnd,rndnorm
real(8),allocatable :: varall(:),stdall(:),thetaPar(:)

z=0
u=rndunif()
if(u.le.cdfPriorKappa(1)) then
  z=1
else if(u.gt.cdfPriorKappa(MaxNCos+1)) then
  z=MaxNCos+1
else
  do i=2,(MaxNCos+1)
    if(cdfPriorKappa(i-1).lt.u .and. u.le.cdfPriorKappa(i)) then
      z=i
      exit
    end if
  end do
end if
kappaPar=KappaGrid(z)

ctol=0.00001d0

gmin=-dlog(ctol)/dble(MaxNCos)
gamPar=gmin+dlog(rndunif()*(dexp(w0*(gmax-gmin))-1.d0)+1.d0)/w0

tau2Par=v0/(2.d0*gamrnd(u0/2.d0,1.d0))
tauPar=dsqrt(tau2Par)
ymean=matmul(dmat,betaPar)
loglike12=sum(matmul(ddata,betaPar))
if(kappaPar.gt.0) then
  allocate(varall(kappaPar),stdall(kappaPar),thetaPar(kappaPar))
  varall=tau2*dexp(-kall(1:kappaPar)*gamPar)
  stdall=dsqrt(varall)
  thetaPar=(/ (rndnorm(), i=1,kappaPar) /)
  thetaPar=stdall*thetaPar
  ymean=ymean+matmul(phi(:,1:kappaPar),thetaPar)
  loglike12=loglike12+sum(matmul(phidata(:,1:kappaPar),thetaPar))
  deallocate(varall,stdall,thetaPar)
end if
c=maxval(ymean)-3.d0
call intsim(dexp(ymean-c),xdelta,nint,fint)
loglike12=loglike12-dble(nobs)*(c+dlog(fint))

infp(1)=ExactLogLikeParOld
infp(2)=loglike12

probz=dexp(infp-maxval(infp)+3.d0)*PriorProbs
probz=probz/sum(probz)
Probzz(1)=probz(1)

ymean=matmul(dmat,beta)
c=maxval(ymean)-3.d0
call intsim(dexp(ymean-c),xdelta,nint,fint)
loglike21=sum(matmul(ddata,beta))-dble(nobs)*(c+dlog(fint))

infp(1)=loglike21
infp(2)=ExactLogLikeSemiOld

probz=dexp(infp-maxval(infp)+3.d0)*PriorProbs
probz=probz/sum(probz)
Probzz(2)=probz(2)

return
end subroutine ProbStay

end subroutine bsad
