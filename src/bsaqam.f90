subroutine bsaqam(verbose,yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,p,&
                  theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,&
                  alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,&
                  sigma2_m0,sigma2_v0,iflagprior,iflagpsi,&
                  maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,&
                  zetag,tau2g,gammag,thetag,betag,alphag,sigma2g,psig,omegag,nug,&
                  fxgridg,fxobsg,yestg,wbg,loglikeg,logpriorg,imodmetg,pmetg)
use ToolsRfunf
use bsamTools
implicit none

!input arguments
integer,intent(in) :: nobs,nparw,nfun,nbasis,nint,fmodel(nfun),iflagprior,iflagpsi
integer,intent(in) :: maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,verbose
real(8), intent(in) :: yobs(nobs),wdata(nobs,nparw),xobs(nobs,nfun),fpm(nfun),p
real(8), intent(in) :: theta0_m0,theta0_s0
real(8), intent(in) :: tau2_m0,tau2_v0,w0,beta_m0(nparw),beta_v0(nparw,nparw)
real(8), intent(in) :: alpha_m0,alpha_s0,psi_m0,psi_s0,omega_m0(nfun),omega_s0(nfun)
real(8), intent(in) :: psifixed,sigma2_m0,sigma2_v0

!output arguments
integer,intent(out) :: imodmetg
real(8),intent(out) :: zetag(smcmc,nfun),tau2g(smcmc,nfun),gammag(smcmc,nfun)
real(8),intent(out) :: thetag(nbasis+1,nfun,smcmc),betag(smcmc,nparw)
real(8),intent(out) :: alphag(smcmc,nfun),sigma2g(smcmc),psig(smcmc,nfun)
real(8),intent(out) :: omegag(smcmc,nfun),fxgridg(nint+1,nfun,smcmc)
real(8),intent(out) :: fxobsg(nobs,nfun,smcmc),wbg(smcmc,nobs),yestg(smcmc,nobs)
real(8),intent(out) :: loglikeg(smcmc),logpriorg(smcmc),pmetg(nfun),nug(smcmc,nobs)

!internal arguments
real(8)  :: stime,itime
integer :: imcmc,isave,nmcmc

integer :: iloop,ifun,iobs
integer :: quadfacts((nbasis+1)*(nbasis+2)/2,3),leftindex((nbasis+1)*(nbasis+2)/2)
integer :: rightindex((nbasis+1)*(nbasis+2)/2),multfact((nbasis+1)*(nbasis+2)/2)
integer :: intsimpfacts(nint+1),tmpA(nbasis+1,1),tmpB(1,nbasis+1)
integer :: tmpI(nbasis+1,nbasis+1),tmpC(nbasis+1,nbasis+1),xinxgrid(nobs,nfun)

real(8) :: xmin(nfun),xmax(nfun),xrange(nfun),xmid(nfun),xdelta(nfun)
real(8) :: wtw(nparw,nparw),wtwi(nparw,nparw),wdatat(nparw,nobs)
real(8) :: xgrid(nint+1,nfun),bhat(nparw),yhat(nobs),shat,xbar(nfun)
real(8) :: xobs2(nobs),xgrid2(nint+1),xidelta(nobs,nfun)

integer :: k,nr
real(8) :: kall(nbasis),kall0(nbasis+1),kbar,kvec(nbasis)

integer :: nfree,nfunconstraint,ifree,irest
real(8),allocatable :: phixobsfree(:,:,:),phixobsfreet(:,:,:)
real(8),allocatable :: phixgridfree(:,:,:),phi2(:,:,:)
real(8),allocatable :: phixobs(:,:,:),phixgrid(:,:,:)

real(8) :: theta0_v0,theta0_v0i,theta0_v0im0,theta0_lnv0
real(8) :: tau2_r0,tau2_s0,tau2_rn,tau2_u0,wk
real(8) :: beta_v0i(nparw,nparw),beta_v0im0(nparw),beta_lnv0
real(8) :: sigma2_r0,sigma2_s0,sigma2_rn,alpha_v0,alpha_v0i,alpha_v0im0,alpha_lnv0
real(8) :: psi_v0i,psi_v0,psi_lnp0,psi_lnv0
real(8) :: omega_v0i(nfun),omega_v0(nfun),omega_lnp0(nfun),omega_lnv0(nfun)

integer :: imodmet,pmet(nfun),iflag_AM,icount_AM,pok
real(8) :: metw,met_alpha,met_omega_alpha,met_psi_alpha
real(8) :: metm(nfun),mets(nfun),metv(nfun),met_beta(nfun),met_beta_AM(nfun)
real(8) :: met_mAM(nfun),met_vAM(nfun),met_omega_m(nfun),met_omega_beta(nfun)
real(8) :: met_psi_m(nfun),met_psi_beta(nfun),met_var_all(nfun)

real(8) :: eta1,eta1sq,eta2sq,nu(nobs),nu_v(nobs,nobs),nu_vi(nobs,nobs)

real(8) :: sigma,sigma2,tau(nfun),tau2(nfun),tau2i(nfun),gampar(nfun),lngampar(nfun)
real(8) :: theta(nbasis+1,nfun),theta02(nfun),theta2(nbasis,nfun),beta(nparw),wb(nobs)
real(8) :: theta0(nfun),gamvec0(nbasis+1,nfun),thv0(nbasis,nfun),thv00(nbasis+1,nfun)
real(8) :: psi(nfun),omega(nfun),zeta(nfun),alpha(nfun),gamvec(nbasis,nfun)
real(8) :: fxobs(nobs,nfun),fxgrid(nint+1,nfun),yest(nobs),ehat(nobs),sigma2i

real(8) :: cdfnorm,rndnorm

call cpu_time(stime)
call rndstart()

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

metw=0.5d0
metm=0.01d0
met_alpha=3.d0
if (met_alpha.gt.2.d0) then
  mets=metm/dsqrt(met_alpha-2.d0)
  metv=mets**2.d0
else
  mets=1.d0
  metv=1.d0
end if
met_beta=(met_alpha-1.d0)*metm
met_beta_AM=met_beta
met_var_all=metm
met_mAM=metm
met_vAM=metv

met_omega_m=0.0001d0
met_omega_alpha=3.d0
met_omega_beta=(met_omega_alpha-1.d0)*met_omega_m

met_psi_m=0.0001d0
met_psi_alpha=3.d0
met_psi_beta=(met_psi_alpha-1.d0)*met_psi_m

wdatat=transpose(wdata)
wtw=matmul(wdatat,wdata)
call inverse(wtw,nparw,wtwi)
bhat=matmul(wtwi,matmul(wdatat,yobs))
yhat=matmul(wdata,bhat)
shat=dsqrt(sum((yobs-yhat)**2.d0)/dble(nobs))
xbar=sum(xobs,1)/dble(nobs)
do ifun=1,nfun
  xmin(ifun)=minval(xobs(:,ifun))
  xmax(ifun)=maxval(xobs(:,ifun))
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

nr=(nbasis+1)*(nbasis+2)/2
kall=(/ (dble(k),k=1,nbasis) /)
kall0=(/ (dble(k),k=0,nbasis) /)
kbar=sum(kall0)/dble(nbasis+1)
kvec=w0/(w0+kall)

eta1=(1.d0-2.d0*p)/(p*(1.d0-p))
eta1sq=eta1**2.d0
eta2sq=2.d0/(p*(1.d0-p))

theta0_v0=theta0_s0**2.d0
theta0_v0i=1.d0/theta0_v0
theta0_v0im0=theta0_v0i*theta0_m0
theta0_lnv0=dlog(theta0_v0)

tau2_r0=2.d0*(2.d0+tau2_m0**2.d0/tau2_v0)
tau2_s0=tau2_m0*(tau2_r0-2.d0)
tau2_rn=tau2_r0+dble(nbasis)

tau2_u0=1/tau2_m0

wk=sum(kall)/2.d0-w0

call inverse(beta_v0,nparw,beta_v0i)
beta_v0im0=matmul(beta_v0i,beta_m0)
beta_lnv0=dlog(determinant(beta_v0,nparw))

alpha_v0=alpha_s0**2.d0
alpha_v0i=1.d0/alpha_v0
alpha_v0im0=alpha_m0*alpha_v0i
alpha_lnv0=dlog(alpha_v0)

psi_v0=psi_s0**2.d0
psi_v0i=1.d0/psi_v0
psi_lnp0=cdfnorm(-psi_m0/psi_s0,0.d0,1.d0,0,1)
psi_lnv0=dlog(psi_v0)

omega_v0=omega_s0**2.d0
omega_v0i=1.d0/omega_v0
omega_lnv0=dlog(omega_v0)
do ifun=1,nfun
  omega_lnp0(ifun)=dlog(cdfnorm((xmax(ifun)-omega_m0(ifun))/(omega_s0(ifun)),0.d0,1.d0,1,0)- &
                        cdfnorm((xmin(ifun)-omega_m0(ifun))/(omega_s0(ifun)),0.d0,1.d0,1,0))
end do

sigma2_r0=2.d0*(2.d0+sigma2_m0**2.d0/sigma2_v0)
sigma2_s0=sigma2_m0*(sigma2_r0-2.d0)
sigma2_rn=sigma2_r0+3.d0*dble(nobs)+dble(nparw)
do ifun=1,nfun
  if(fmodel(ifun).eq.1) then
    sigma2_rn=sigma2_rn+dble(nbasis)
  else if (fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
    sigma2_rn=sigma2_rn+dble(nbasis+1)/2.d0+1.d0
  else
    sigma2_rn=sigma2_rn+dble(nbasis+1)/2.d0
  end if
end do

sigma=shat
sigma2=sigma**2.d0
sigma2i=1.d0/sigma2
tau2=tau2_m0
tau=dsqrt(tau2)
tau2i=1.d0/tau2
psi=psifixed
omega=xmid
gampar=1.d0/w0
lngampar=dlog(gampar)
do ifun=1,nfun
  gamvec(:,ifun)=dexp(-gampar(ifun)*kall)
  gamvec0(1,ifun)=1.d0
  gamvec0(2:(nbasis+1),ifun)=gamvec(:,ifun)
  thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun)
  thv00(1,ifun)=theta0_v0
  thv00(2:(nbasis+1),ifun)=thv0(:,ifun)
  zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
  theta(1,ifun)=0.1d0
  theta0(ifun)=theta(1,ifun)
  theta02(ifun)=theta0(ifun)**2.d0
  do k=2,nbasis+1
    theta(k,ifun)=0.1d0*dsqrt(gamvec(k-1,ifun))*rndnorm()
    theta2(k-1,ifun)=theta(k,ifun)**2.d0
  end do
end do
beta=0.d0
wb=matmul(wdata,beta)
alpha=0.d0
call diag(1.d0,nobs,nu_v)
call diag(1.d0,nobs,nu_vi)
nu=sigma2

nfree=count(fmodel.eq.1)
nfunconstraint=count(fmodel.gt.1)
allocate(phixobsfree(nobs,nbasis,nfree),phixobs(nr,nobs,nfunconstraint))
allocate(phixgridfree(nint+1,nbasis,nfree),phixgrid(nr,nint+1,nfunconstraint))
allocate(phi2(nbasis,nbasis,nfree),phixobsfreet(nbasis,nobs,nfree))
ifree=1
irest=1
do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobsfree(:,:,ifree))
    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgridfree(:,:,ifree))
    phixobsfreet(:,:,ifree)=transpose(phixobsfree(:,:,ifree))
    phi2(:,:,ifree)=matmul(phixobsfreet(:,:,ifree),phixobsfree(:,:,ifree))

    call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
                  nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))

    ifree=ifree+1
  else
    if (fmodel(ifun).eq.2) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
                  IntConst2,IntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
                  IntConst2,IntCos,nbasis,nint+1,phixgrid(:,:,irest))

      call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest), &
                  quadfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))

    end if

    if (fmodel(ifun).eq.3) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))

      call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                      xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                      nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    end if

    if (fmodel(ifun).eq.4) then
      xobs2=xmin(ifun)+xmax(ifun)-xobs(:,ifun)
      xgrid2=xmin(ifun)+xmax(ifun)-xgrid(:,ifun)
      call GetPhi(xobs2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
                  IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))

      call GetConcavef(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                       xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                       nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    end if

    if (fmodel(ifun).eq.5) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))

      call GetSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),&
                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    end if

    if (fmodel(ifun).eq.6) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))

      call GetRotateSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun), &
                       xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
                       xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
                       xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
                       fxobs(:,ifun),fxgrid(:,ifun))
    end if

    if (fmodel(ifun).eq.7) then
      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))

      call GetUf(fpm(ifun),omega(ifun),psi(ifun),theta(:,ifun),xobs(:,ifun), &
                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),quadfacts, &
                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    end if

    irest=irest+1
  end if
end do

if(maxval(fmodel).gt.1) then
  if (verbose.eq.1) then
    call dblepr('Initializing MCMC parameters ...',-1,1.d0,0)
  end if
  do imodmet=1,maxmodmet
    imodmetg=imodmet

    pmet=0
    iflag_AM=0
    icount_AM=0
    do imcmc=1,nblow0
      call rchkusr()
      call GetMCMC()
    end do

    pmet=0
    iflag_AM=1
    icount_AM=0
    do imcmc=1,nblow0
      call rchkusr()
      call GetMCMC()
    end do

    pok=0
    do ifun=1,nfun
      if (fmodel(ifun).gt.1) then
        if (dble(pmet(ifun))/dble(nblow0).gt.0.6d0) then
          if (verbose.eq.1) then
            call haprint(ifun, dble(pmet(ifun))/dble(nblow0))
          end if

          metm(ifun)=metm(ifun)*10.d0
          met_var_all(ifun)=metm(ifun)
          met_mAM(ifun)=metm(ifun)
          if (met_alpha.gt.2.d0) then
            mets(ifun)=metm(ifun)/dsqrt(met_alpha-2.d0)
            metv(ifun)=mets(ifun)**2.d0
          else
            mets(ifun)=1.d0
            metv(ifun)=1.d0
          end if
          met_vAM(ifun)=metv(ifun)

          met_omega_m(ifun)=10.d0*met_omega_m(ifun)
          met_omega_beta(ifun)=(met_omega_alpha-1.d0)*met_omega_m(ifun)

          met_psi_m(ifun)=10.d0*met_psi_m(ifun)
          met_psi_beta(ifun)=(met_psi_alpha-1.d0)*met_psi_m(ifun)

          theta(1,ifun)=0.1d0
          do k=2,(nbasis+1)
            theta(k,ifun)=0.01d0*dsqrt(sigma)*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
          end do
          tau2(ifun)=tau2_m0
          tau(ifun)=dsqrt(tau2_m0)
          gampar(ifun)=1.d0/w0
          lngampar(ifun)=dlog(gampar(ifun))
          gamvec(:,ifun)=dexp(-kall*gampar(ifun))
          psi(ifun)=psifixed
          omega(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
          zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
          alpha(ifun)=0.d0
        else if (dble(pmet(ifun))/dble(nblow0).lt.0.3d0) then
          if (verbose.eq.1) then
            call laprint(ifun, dble(pmet(ifun))/dble(nblow0))
          end if

          metm(ifun)=metm(ifun)/10.d0
          met_var_all(ifun)=metm(ifun)
          met_mAM(ifun)=metm(ifun)
          if (met_alpha.gt.2.d0) then
            mets(ifun)=metm(ifun)/dsqrt(met_alpha-2.d0)
            metv(ifun)=mets(ifun)**2.d0
          else
            mets(ifun)=1.d0
            metv(ifun)=1.d0
          end if
          met_vAM(ifun)=metv(ifun)

          met_omega_m(ifun)=met_omega_m(ifun)/10.d0
          met_omega_beta(ifun)=(met_omega_alpha-1.d0)*met_omega_m(ifun)

          met_psi_m(ifun)=met_psi_m(ifun)/10.d0
          met_psi_beta(ifun)=(met_psi_alpha-1.d0)*met_psi_m(ifun)

          theta(1,ifun)=0.1d0
          do k=2,(nbasis+1)
            theta(k,ifun)=0.01d0*dsqrt(sigma)*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
          end do
          tau2(ifun)=tau2_m0
          tau(ifun)=dsqrt(tau2_m0)
          gampar(ifun)=1.d0/w0
          lngampar(ifun)=dlog(gampar(ifun))
          gamvec(:,ifun)=dexp(-kall*gampar(ifun))
          psi(ifun)=psifixed
          omega(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
          zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
          alpha(ifun)=0.d0
        else
          pok=pok+1
        end if
      end if
    end do

    if (pok.eq.nfunconstraint) then
      exit
    end if

    if (imodmet.lt.maxmodmet) then
      if (met_alpha.gt.2.d0) then
        mets=metm/dsqrt(met_alpha-2.d0)
        metv=mets**2.d0
      else
        mets=1.d0
        metv=1.d0
      end if
      met_beta=(met_alpha-1.d0)*metm
      met_beta_AM=met_beta
      met_var_all=metm
      met_mAM=metm
      met_vAM=metv

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
    end if
  end do
end if

isave=1
nmcmc=nblow+nskip*smcmc
do imcmc=1,nmcmc
  if(imcmc.eq.1) then
    if (verbose.eq.1) then
      call dblepr('Burnin ...',-1,1.d0,0)
    end if
    pmet=0
  end if
  call rchkusr()
  call GetMCMC()

  if(imcmc.eq.nblow) then
    if (verbose.eq.1) then
      do ifun=1,nfun
        if(fmodel(ifun).gt.1) then
          call aprint(ifun, dble(pmet(ifun))/dble(nblow))
        end if
      end do
      call dblepr('Main iterations ...',-1,1.d0,0)
    end if
    pmet=0
  end if

  if(imcmc.gt.nblow .and. mod(imcmc,nskip).eq.0) then
    betag(isave,:)=beta
    sigma2g(isave)=sigma2
    zetag(isave,:)=zeta
    tau2g(isave,:)=tau2
    gammag(isave,:)=gampar
    thetag(:,:,isave)=theta
    alphag(isave,:)=alpha
    psig(isave,:)=psi
    omegag(isave,:)=omega
    nug(isave,:)=nu

    wbg(isave,:)=wb
    fxobsg(:,:,isave)=fxobs
    fxgridg(:,:,isave)=fxgrid
    yest=wb+sum(fxobs,2)
    yestg(isave,:)=yest

    ehat=yobs-yest
    loglikeg(isave)=-sum(dabs(ehat)+(2.d0*p-1.d0)*ehat)/(2.d0*sigma2)+ &
                    dble(nobs)*(dlog(p)+dlog(1.d0-p)-dlog(sigma2))

    logpriorg(isave)=GetLogPrior()

    if (verbose.eq.1) then
      if (mod(isave,ndisp).eq.0) then
        call cpu_time(itime)
        call sprint(isave,smcmc,itime-stime)
      end if
    end if
    isave=isave+1
  end if
end do
deallocate(phixobsfree,phixobsfreet,phixobs,phixgridfree,phixgrid,phi2)
pmetg=dble(pmet)/dble(smcmc*nskip)
if (verbose.eq.1) then
  do ifun=1,nfun
    if(fmodel(ifun).gt.1) then
      call aprint(ifun, pmetg(ifun))
    end if
  end do
end if
call rndend()

!=========================================================================================
contains

subroutine GetMCMC()
implicit none

!Internal arguments
real(8) :: gamrnd,rndunif,rtgamrnd,normrnd,ltnormrnd
real(8) :: rtnormrnd,tnormrnd,ltgamrnd,invgaussrnd
real(8) :: resid(nobs),vni(nparw,nparw),vn(nparw,nparw),bn(nparw)
real(8) :: lambda1(nobs),lambda2
real(8) :: residb(nparw),sse,ck,c2
real(8) :: rk(nobs),vni1(nbasis,nbasis),vn1(nbasis,nbasis),bn1(nbasis)
real(8) :: met_var_all_new(nfun),met_beta_new(nfun),testp
real(8) :: theta0_old,theta_old(nbasis),met_var0,met_var,met_std0
real(8) :: theta0_new,theta02_new,theta_new(nbasis),theta2_new(nbasis),thetanew(nbasis+1)
real(8) :: theta0_new_lnp,theta0_old_lnp,fxobs_new(nobs),fxgrid_new(nint+1)
real(8) :: met_stdS,met_varS,psi_old,psi_lnpold,psi_new,psi_lnpnew
real(8) :: omega_old,omega_lnpold,omega_new,omega_lnpnew
real(8) :: resid_old(nobs),resid_new(nobs),sse_old,sse_new,sold,snew
real(8) :: fx1(nobs),fxg1(nint+1),a_vni,a_vn,a_mn,xtx(nfun)
real(8) :: th2gam(nbasis),sn,bup
real(8) :: ck1(nbasis),ckz(nbasis),u1(nbasis),u2,bnz(nbasis),bmin
real(8),allocatable :: bnz1(:)
integer :: z(nbasis)

resid=yobs-sum(fxobs,2)-eta1*nu
vni=matmul(wdatat,matmul(nu_vi,wdata))/eta2sq+beta_v0i
call inverse(vni,nparw,vn)
bn=beta_v0im0+matmul(wdatat,matmul(nu_vi,resid))/eta2sq
bn=matmul(vn,bn)

call mvnrnd(bn,sigma2*vn,nparw,beta)
wb=matmul(wdata,beta)

resid=yobs-wb-sum(fxobs,2)
lambda2=2.d0/sigma2+eta1sq/(sigma2*eta2sq)
lambda1=dsqrt((eta1sq+2.d0*eta2sq)/(resid**2.d0))
do iobs=1,nobs
  nu(iobs)=1.d0/invgaussrnd(lambda1(iobs),lambda2)
end do
call diagvec(nu,nobs,nu_v)
call diagvec(1.d0/nu,nobs,nu_vi)

resid=yobs-wb-sum(fxobs,2)-eta1*nu
sse=sigma2_s0+sum(resid**2.d0/nu)/eta2sq

sse=sse+2.d0*sum(nu)

residb=beta-beta_m0
sse=sse+dot_product(residb,matmul(beta_v0i,residb))

if (maxval(fmodel).eq.1) then
  sse=sse+sum(theta2/thv0)
  sigma2=sse/(2.d0*gamrnd(sigma2_rn/2.d0,1.d0))
  sigma=dsqrt(sigma2)
  sigma2i=1.d0/sigma2
else
  ck=0.d0
  do ifun=1,nfun
    if (fmodel(ifun).eq.1) then
      sse=sse+sum(theta2(:,ifun)/thv0(:,ifun))
    else
      ck=ck+theta02(ifun)/(2.d0*theta0_v0)+sum(theta2(:,ifun)/thv0(:,ifun))/2.d0
    end if

    if (fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
      sse=sse+((alpha(ifun)-alpha_m0)**2.d0)/alpha_v0
    end if
  end do
  c2=1.d0/sigma-dlog(rndunif())/ck
  c2=c2**2.d0
  sigma2i=rtgamrnd(sigma2_rn/2.d0,2.d0/sse,c2)
  sigma2=1.d0/sigma2i
  sigma=dsqrt(sigma2)
end if

do ifun=1,nfun
  gamvec0(1,ifun)=1.d0
  gamvec0(2:(nbasis+1),ifun)=gamvec(:,ifun)
  thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun)
  thv00(1,ifun)=theta0_v0
  thv00(2:(nbasis+1),ifun)=thv0(:,ifun)
  theta0(ifun)=theta(1,ifun)
  theta02(ifun)=theta0(ifun)**2.d0
  theta2(:,ifun)=theta(2:(nbasis+1),ifun)**2.d0
end do

met_var_all_new=0.d0
met_beta_new=0.d0
resid=yobs-wb-sum(fxobs,2)-eta1*nu
ifree=1
irest=1
do ifun=1,nfun
  testp=0.d0
  rk=resid+fxobs(:,ifun)
  if (fmodel(ifun).eq.1) then
    vni1=matmul(phixobsfreet(:,:,ifree),matmul(nu_vi,phixobsfree(:,:,ifree)))/eta2sq
    do k=1,nbasis
      vni1(k,k)=vni1(k,k)+1.d0/thv0(k,ifun)
    end do
    call inverse(vni1,nbasis,vn1)
    bn1=matmul(vn1,matmul(phixobsfreet(:,:,ifree),matmul(nu_vi,rk))/eta2sq)
    call mvnrnd(bn1,sigma2*vn1,nbasis,theta(2:(nbasis+1),ifun))
    theta(1,ifun)=0.d0
    theta02(ifun)=0.d0
    theta2(:,ifun)=theta(2:(nbasis+1),ifun)**2.d0
    call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
                  nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))

    ifree=ifree+1
  else
    theta0_old=theta(1,ifun)
    theta_old=theta(2:(nbasis+1),ifun)

    met_var_all_new(ifun)=met_beta_AM(ifun)/gamrnd(met_alpha,1.d0)
    ck=metw*met_var_all_new(ifun) + (1.d0-metw)*met_mAM(ifun)
    met_beta_new(ifun)=(met_alpha-1.d0)*ck

    ck=5.66d0
    met_var0=ck*met_var_all_new(ifun)
    met_var=ck*met_var_all_new(ifun)

    met_std0=dsqrt(met_var0*sigma)
    theta0_new=ltnormrnd(theta0_old,met_std0,0.d0)
    do k=1,nbasis
      theta_new(k)=normrnd(theta_old(k),dsqrt(met_var*thv0(k,ifun)*sigma))
    end do

    theta02_new=theta0_new**2.d0
    theta2_new=theta_new**2.d0

    theta0_new_lnp=cdfnorm(-theta0_old/met_std0,0.d0,1.d0,0,1)
    theta0_old_lnp=cdfnorm(-theta0_new/met_std0,0.d0,1.d0,0,1)

    thetanew(1)=theta0_new
    thetanew(2:(nbasis+1))=theta_new
    if (fmodel(ifun).eq.2) then
      call  GetUpf(fpm(ifun),thetanew,phixobs(:,:,irest),phixgrid(:,:,irest),&
                   quadfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
    else if (fmodel(ifun).eq.3) then
      call GetConvexf(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
                      xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                      nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
    else if (fmodel(ifun).eq.4) then
      call GetConcavef(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
                       xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
                       nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
    else
      if (iflagpsi.eq.1) then
        met_varS=met_psi_beta(ifun)/gamrnd(met_psi_alpha,1.d0)
        met_stdS=dsqrt(met_varS)

        psi_old=psi(ifun)
        psi_new=ltnormrnd(psi_old,met_stdS,0.d0)
        psi_lnpnew=cdfnorm(-psi_old/met_stdS,0.d0,1.d0,0,1)
        psi_lnpold=cdfnorm(-psi_new/met_stdS,0.d0,1.d0,0,1)

        testp=testp-&
              ((psi_new-psi_m0)**2.d0)/(2.d0*psi_v0)+&
              ((psi_old-psi_m0)**2.d0)/(2.d0*psi_v0)-&
              psi_lnpold+ &
              psi_lnpnew
      else
        psi_old=psi(ifun)
        psi_new=psi(ifun)
      end if

      met_varS=met_omega_beta(ifun)/gamrnd(met_omega_alpha,1.d0)
      met_stdS=dsqrt(met_varS)

      omega_old=omega(ifun)
      omega_new=tnormrnd(omega_old,met_stdS,xmin(ifun),xmax(ifun))
      omega_lnpnew=dlog(cdfnorm((xmax(ifun)-omega_old)/met_stdS,0.d0,1.d0,1,0)- &
                        cdfnorm((xmin(ifun)-omega_old)/met_stdS,0.d0,1.d0,1,0))
      omega_lnpold=dlog(cdfnorm((xmax(ifun)-omega_new)/met_stdS,0.d0,1.d0,1,0)- &
                        cdfnorm((xmin(ifun)-omega_new)/met_stdS,0.d0,1.d0,1,0))

      testp=testp- &
            ((omega_new-omega_m0(ifun))**2.d0)/(2.d0*omega_v0(ifun))+&
            ((omega_old-omega_m0(ifun))**2.d0)/(2.d0*omega_v0(ifun))-&
            omega_lnpold+ &
            omega_lnpnew

      if (fmodel(ifun).eq.5) then
        call GetSf(fpm(ifun),omega_new,psi_new,alpha(ifun),thetanew,xobs(:,ifun),&
                   xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
                   xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
                   intsimpfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
      else if (fmodel(ifun).eq.6) then
        call GetRotateSf(fpm(ifun),omega_new,psi_new,alpha(ifun),thetanew,xobs(:,ifun),&
                         xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
                         xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
                         xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
                         fxobs_new,fxgrid_new)
      else
        call GetUf(fpm(ifun),omega_new,psi_new,thetanew,xobs(:,ifun), &
                   xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
                   xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),quadfacts, &
                   intsimpfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
      end if
    end if

    resid_new=rk-fxobs_new
    sse_new=sum(resid_new**2.d0/nu)

    resid_old=rk-fxobs(:,ifun)
    sse_old=sum(resid_old**2.d0/nu)

    snew=sum(theta2_new/thv0(:,ifun))
    sold=sum(theta2(:,ifun)/thv0(:,ifun))

    testp=testp- &
          (sse_new-sse_old)/(2.d0*sigma2*eta2sq)- &
          (theta02_new)/(2.d0*sigma*theta0_v0)+ &
          (theta02(ifun))/(2.d0*sigma*theta0_v0)- &
          (snew-sold)/(2.d0*sigma)- &
          theta0_old_lnp+theta0_new_lnp- &
          ((met_alpha+1.d0)*dlog(met_var_all(ifun)))-&
          (met_beta_new(ifun)/met_var_all(ifun))+ &
          ((met_alpha+1.d0)*dlog(met_var_all_new(ifun)))+&
          (met_beta_AM(ifun)/met_var_all_new(ifun))

    if(dlog(rndunif()).le.testp) then
      theta0(ifun)=theta0_new
      theta(1,ifun)=theta0_new
      theta(2:(nbasis+1),ifun)=theta_new
      theta2(:,ifun)=theta2_new
      theta02(ifun)=theta02_new
      if (fmodel(ifun).eq.5 .or. fmodel(ifun).eq.6 .or. fmodel(ifun).eq.7) then
        psi(ifun)=psi_new
        omega(ifun)=omega_new
      end if
      met_var_all(ifun)=met_var_all_new(ifun)
      fxobs(:,ifun)=fxobs_new
      fxgrid(:,ifun)=fxgrid_new
      pmet(ifun)=pmet(ifun)+1
    end if

    if (iflag_AM.eq.1) then
      icount_AM=icount_AM+1
      met_mAM(ifun)=met_mAM(ifun)+(met_var_all(ifun)-met_mAM(ifun))/dble(icount_AM)
      met_vAM(ifun)=(dble(icount_AM-1)/dble(icount_AM))*met_vAM(ifun) + &
                    ((met_var_all(ifun)-met_mAM(ifun))**2.d0)/dble(icount_AM)
    end if
    ck=metw*met_var_all(ifun)+(1.d0-metw)*met_mAM(ifun)
    met_beta_AM(ifun)=(met_alpha-1.d0)*ck

    irest=irest+1
  end if
  resid=yobs-wb-sum(fxobs,2)-eta1*nu
end do

do ifun=1,nfun
  if(fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
    resid=yobs-wb-sum(fxobs,2)-eta1*nu
    fx1=fxobs(:,ifun)-alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxg1=fxgrid(:,ifun)-alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))

    rk=resid+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))

    xtx(ifun)=sum((xobs(:,ifun)-xmid(ifun))**2.d0/nu)
    a_vni=xtx(ifun)/eta2sq+alpha_v0i
    a_vn=1.d0/a_vni
    a_mn=a_vn*(sum((xobs(:,ifun)-xmid(ifun))*rk/nu)/eta2sq+alpha_v0im0)
    if (fpm(ifun).eq.1.d0) then
      alpha(ifun)=ltnormrnd(a_mn,dsqrt(sigma2*a_vn),0.d0)
    else
      alpha(ifun)=rtnormrnd(a_mn,dsqrt(sigma2*a_vn),0.d0)
    end if
    fxobs(:,ifun)=fx1+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxgrid(:,ifun)=fxg1+alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))
  else
    alpha(ifun)=0.d0
  end if
end do

do ifun=1,nfun
  th2gam=theta2(:,ifun)/gamvec(:,ifun)
  if (fmodel(ifun).eq.1) then
    th2gam=th2gam/sigma2
  else
    th2gam=th2gam/sigma
  end if
  sn=sum(th2gam)
  if (iflagprior.eq.0) then
    tau2(ifun)=(tau2_s0+sn)/(2.d0*gamrnd(tau2_rn/2.d0,1.d0))
    tau(ifun)=dsqrt(tau2(ifun))
  else
    bup=tau2(ifun)-dlog(rndunif())/tau2_u0

    tau2i(ifun)=ltgamrnd(dble(nbasis-2)/2.d0,2.d0/sn,1.d0/bup);
    tau2(ifun)=1.d0/tau2i(ifun)
    tau(ifun)=dsqrt(tau2(ifun))
  end if
  thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun)

  if (fmodel(ifun).eq.1) then
    ck1=theta2(:,ifun)/(2.d0*sigma2*tau2(ifun))
  else
    ck1=theta2(:,ifun)/(2.d0*sigma*tau2(ifun))
  end if

  if (count(ck1.eq.0.d0).eq.nbasis) then
    gampar(ifun) = 1.d0/w0
  else
    z=transfer(ck1.eq.0.d0, (/ (1,k=1,nbasis) /))
    ckz=dble(z)*(1.d0-ck1)+ck1
    u1=(/ (rndunif(),k=1,nbasis) /)
    bnz=gampar(ifun)+(dlog(ckz-dlog(u1)*gamvec(:,ifun))-dlog(ckz))/kall
    if (count(z.eq.1).gt.0) then
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
    gampar(ifun)=bmin+dlog(u2+(1.d0-u2)*dexp(-wk*bmin))/wk
  end if

  gamvec(:,ifun)=dexp(-gampar(ifun)*kall)
  lngampar(ifun)=dlog(gampar(ifun))
  zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
end do

return
end subroutine GetMCMC


function GetLogPrior()
implicit none

!Output argument
real(8) :: GetLogPrior

!Internal argument
integer :: ifun
real(8) :: residb(nparw),thvar(nbasis+1),thetasq(nbasis+1,nfun)
real(8) :: gamvec(nbasis,nfun),sse,theta0_lnp0,alpha_lnp0,lnpriorf
real(8) :: cdfnorm

lnpriorf=0.d0
thetasq=theta**2.d0
do ifun=1,nfun
  gamvec(:,ifun)=dexp(-kall*gampar(ifun))
end do

residb=beta-beta_m0
lnpriorf=lnpriorf-dot_product(residb,matmul(beta_v0i,residb))/(2.d0*sigma2)- &
         dble(nparw)*dlog(2.d0*PI*sigma2)/2.d0-beta_lnv0/2.d0

lnpriorf=lnpriorf+LogfIG(sigma2,sigma2_r0,sigma2_s0)

do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    thvar(2:(nbasis+1))=sigma2*tau2(ifun)*gamvec(:,ifun)
    sse=sum(thetasq(2:(nbasis+1),ifun)/thvar(2:(nbasis+1)))
    lnpriorf=lnpriorf-sse/2.d0-dble(nbasis)*dlog(2.d0*PI)/2.d0- &
             sum(dlog(thvar(2:(nbasis+1))))/2.d0

  else
    thvar(1)=sigma*theta0_v0
    thvar(2:(nbasis+1))=sigma*tau2(ifun)*gamvec(:,ifun)
    sse=sum(thetasq(:,ifun)/thvar)
    theta0_lnp0=cdfnorm(-theta0_m0/(theta0_s0*dsqrt(sigma)),0.d0,1.d0,0,1)
    lnpriorf=lnpriorf-sse/2.d0-dble(nbasis+1)*dlog(2.d0*PI)/2.d0- &
             sum(dlog(thvar))/2.d0-theta0_lnp0
  end if

  lnpriorf=lnpriorf-w0*gampar(ifun)-dlog(w0)

  if (iflagprior.eq.1) then
    lnpriorf=lnpriorf-tau2_u0*tau2(ifun)-tau2_u0
  else
    lnpriorf=lnpriorf+LogfIG(tau2(ifun),tau2_r0,tau2_s0)
  end if

  if (fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
    alpha_lnp0=cdfnorm(-alpha_m0/(sigma*alpha_s0),0.d0,1.d0,0,1)
    lnpriorf=lnpriorf-((alpha(ifun)-alpha_m0)**2.d0)/(2.d0*sigma2*alpha_v0)- &
             dlog(2.d0*PI*alpha_v0)/2.d0-alpha_lnp0
  end if

  if (fmodel(ifun).eq.5 .or. fmodel(ifun).eq.6 .or. fmodel(ifun).eq.7) then
    if (iflagpsi.eq.1) then
      lnpriorf=lnpriorf-((psi(ifun)-psi_m0)**2.d0)/(2.d0*psi_v0)- &
               dlog(2.d0*PI*psi_v0)/2.d0-psi_lnp0
    end if
    lnpriorf=lnpriorf-((omega(ifun)-omega_m0(ifun))**2.d0)/(2.d0*omega_v0(ifun))- &
             dlog(2.d0*PI*omega_v0(ifun))/2.d0-omega_lnp0(ifun)
  end if
end do

GetLogPrior=lnpriorf

return
end function GetLogPrior


end subroutine bsaqam
