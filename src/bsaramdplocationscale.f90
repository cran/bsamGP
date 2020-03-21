subroutine bsaramdplocationscale(verbose,yobs,wdata,xobs,nobs,nparw,ndimw,nfun,nbasis,nint,fmodel,fpm,&
                                 theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,&
                                 alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,&
                                 kappa_r0,kappa_s0,sigma2_r0,sigma2_s0,&
                                 tmass_a,tmass_b,egrid,ngrid,iflagprior,iflagpsi,&
                                 maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,&
                                 zetag,tau2g,gammag,thetag,betag,alphag,psig,omegag,&
                                 fxgridg,fxobsg,yestg,wbg,tmassg,config,nclassg,mug,&
                                 mu0g,kappag,sigma2g,edensg,invlikeg,imodmetg,pmetg)
use ToolsRfunf
use bsamTools
implicit none

!input arguments
integer,intent(in) :: nobs,nparw,ndimw,nfun,nbasis,nint,fmodel(nfun)
integer,intent(in) :: maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,verbose
integer,intent(in) :: ngrid,iflagprior,iflagpsi
real(8), intent(in) :: yobs(nobs),wdata(nobs,nparw),xobs(nobs,nfun),fpm(nfun)
real(8), intent(in) :: theta0_m0,theta0_s0,kappa_r0,kappa_s0,egrid(ngrid)
real(8), intent(in) :: tau2_m0,tau2_v0,w0,beta_m0(nparw),beta_v0(nparw,nparw)
real(8), intent(in) :: alpha_m0,alpha_s0,psi_m0,psi_s0,omega_m0(nfun),omega_s0(nfun)
real(8), intent(in) :: psifixed,sigma2_r0,sigma2_s0,tmass_a,tmass_b

!output arguments
integer,intent(out) :: imodmetg,nclassg(smcmc),config(smcmc,nobs)
real(8),intent(out) :: zetag(smcmc,nfun),tau2g(smcmc,nfun),gammag(smcmc,nfun)
real(8),intent(out) :: thetag(nbasis+1,nfun,smcmc),betag(smcmc,nparw)
real(8),intent(out) :: alphag(smcmc,nfun),sigma2g(smcmc,nobs),psig(smcmc,nfun)
real(8),intent(out) :: omegag(smcmc,nfun),fxgridg(nint+1,nfun,smcmc)
real(8),intent(out) :: fxobsg(nobs,nfun,smcmc),wbg(smcmc,nobs),yestg(smcmc,nobs)
real(8),intent(out) :: invlikeg(smcmc,nobs),pmetg(nfun),tmassg(smcmc)
real(8),intent(out) :: mug(smcmc,nobs),mu0g(smcmc),kappag(smcmc),edensg(smcmc,ngrid)

!internal arguments
real(8)  :: stime,itime
integer :: iadapt,imcmc,isave,nmcmc

integer :: iloop,ifun,iobs
integer :: quadfacts((nbasis+1)*(nbasis+2)/2,3),leftindex((nbasis+1)*(nbasis+2)/2)
integer :: rightindex((nbasis+1)*(nbasis+2)/2),multfact((nbasis+1)*(nbasis+2)/2)
integer :: intsimpfacts(nint+1),tmpA(nbasis+1,1),tmpB(1,nbasis+1)
integer :: tmpI(nbasis+1,nbasis+1),tmpC(nbasis+1,nbasis+1),xinxgrid(nobs,nfun)
real(8) :: xidelta(nobs,nfun)

real(8) :: xmin(nfun),xmax(nfun),xrange(nfun),xmid(nfun),xdelta(nfun)
real(8) :: wtw(nparw,nparw),wtwi(nparw,nparw),wdatat(nparw,nobs)
real(8) :: xgrid(nint+1,nfun),bhat(nparw),yhat(nobs)
real(8) :: xobs2(nobs),xgrid2(nint+1)

integer :: k,nr
real(8) :: kall(nbasis),kall0(nbasis+1),kbar

integer :: nfree,nfunconstraint,ifree,irest
real(8),allocatable :: phixobsfree(:,:,:),phixobsfreet(:,:,:)
real(8),allocatable :: phixgridfree(:,:,:),phixobs(:,:,:),phixgrid(:,:,:)

real(8) :: theta0_v0,theta0_v0i,theta0_v0im0
real(8) :: tau2_r0,tau2_s0,tau2_rn,tau2_u0,wk
real(8) :: beta_v0i(nparw,nparw),beta_v0im0(nparw)
real(8) :: alpha_v0,alpha_v0i,alpha_v0im0
real(8) :: psi_v0i,psi_v0,omega_v0i(nfun),omega_v0(nfun)

integer :: imodmet,pmet(nfun),iflag_AM,icount_AM,pok
real(8) :: metw,met_alpha,met_omega_alpha,met_psi_alpha
real(8) :: metm(nfun),mets(nfun),metv(nfun),met_beta(nfun),met_beta_AM(nfun)
real(8) :: met_mAM(nfun),met_vAM(nfun),met_omega_m(nfun),met_omega_beta(nfun)
real(8) :: met_psi_m(nfun),met_psi_beta(nfun),met_var_all(nfun)

real(8) :: tau(nfun),tau2(nfun),tau2i(nfun),gampar(nfun),lngampar(nfun)
real(8) :: theta(nbasis+1,nfun),theta02(nfun),theta2(nbasis,nfun),beta(nparw),wb(nobs)
real(8) :: theta0(nfun),gamvec0(nbasis+1,nfun),thv0(nbasis,nfun),thv00(nbasis+1,nfun)
real(8) :: psi(nfun),omega(nfun),zeta(nfun),alpha(nfun),gamvec(nbasis,nfun)
real(8) :: fxobs(nobs,nfun),fxgrid(nint+1,nfun),yest(nobs)

integer :: iclass,nclass,classind(nobs)
integer,allocatable :: nclassh(:),e_index(:)
real(8) :: invlike(nobs),sigma2(nobs),sigmai(nobs,nobs)
real(8) :: tmass,edens(ngrid),mu(nobs),mu0,mu0m,kappa
real(8), allocatable :: sigma2h(:),muh(:)

real(8) :: rndnorm,gamrnd

call cpu_time(stime)
call rndstart()

wdatat=transpose(wdata)
wtw=matmul(wdatat,wdata)
call inverse(wtw,nparw,wtwi)
bhat=matmul(wtwi,matmul(wdatat,yobs))
yhat=matmul(wdata,bhat)
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

theta0_v0=theta0_s0**2.d0
theta0_v0i=1.d0/theta0_v0
theta0_v0im0=theta0_v0i*theta0_m0

tau2_r0=2.d0*(2.d0+tau2_m0**2.d0/tau2_v0)
tau2_s0=tau2_m0*(tau2_r0-2.d0)
tau2_rn=tau2_r0+dble(nbasis)

tau2_u0=1/tau2_m0

wk=sum(kall)/2.d0-w0

call inverse(beta_v0,nparw,beta_v0i)
beta_v0im0=matmul(beta_v0i,beta_m0)

alpha_v0=alpha_s0**2.d0
alpha_v0i=1.d0/alpha_v0
alpha_v0im0=alpha_m0*alpha_v0i

psi_v0=psi_s0**2.d0
psi_v0i=1.d0/psi_v0

omega_v0=omega_s0**2.d0
omega_v0i=1.d0/omega_v0

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

nclass=nobs
classind=(/ (iobs,iobs=1,nobs) /)
tmass=1.d0
mu=(/ (rndnorm(),iobs=1,nobs) /)
mu0=0.d0
kappa=1.d0
sigma2=(/ (1.d0/gamrnd(3.d0,1.d0),iobs=1,nobs) /)
call diag(1.d0,nobs,sigmai)

nfree=count(fmodel.eq.1)
nfunconstraint=count(fmodel.gt.1)
allocate(phixobsfree(nobs,nbasis,nfree),phixobs(nr,nobs,nfunconstraint))
allocate(phixgridfree(nint+1,nbasis,nfree),phixgrid(nr,nint+1,nfunconstraint))
allocate(phixobsfreet(nbasis,nobs,nfree))
ifree=1
irest=1
do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobsfree(:,:,ifree))
    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgridfree(:,:,ifree))
    phixobsfreet(:,:,ifree)=transpose(phixobsfree(:,:,ifree))

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

imcmc=0
if(maxval(fmodel).gt.1) then
  if (verbose.eq.1) then
    call piprint()
  end if
  do imodmet=1,maxmodmet
    imodmetg=imodmet

    pmet=0
    iflag_AM=0
    icount_AM=0
    do iadapt=1,nblow0
      call rchkusr()  ! user interrupt
      call GetMCMC()
    end do

    pmet=0
    iflag_AM=1
    icount_AM=0
    do iadapt=1,nblow0
      call rchkusr()  ! user interrupt
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
            theta(k,ifun)=0.01d0*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
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
            theta(k,ifun)=0.01d0*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
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
iflag_AM=1  ! placeholder for Free only case.
nmcmc=nblow+nskip*smcmc
do imcmc=1,nmcmc
  if(imcmc.eq.1) then
    if (verbose.eq.1) then
      call biprint()
    end if
    pmet=0
  end if
  call rchkusr()  ! user interrupt
  call GetMCMC()  ! Draw samples

  if(imcmc.eq.nblow) then
    if (verbose.eq.1) then
      do ifun=1,nfun
        if(fmodel(ifun).gt.1) then
          call aprint(ifun, dble(pmet(ifun))/dble(nblow))
        end if
      end do
      call miprint
    end if
    pmet=0
  end if

  if(imcmc.gt.nblow .and. mod(imcmc,nskip).eq.0) then
    betag(isave,:)=beta
    zetag(isave,:)=zeta
    tau2g(isave,:)=tau2
    gammag(isave,:)=gampar
    thetag(:,:,isave)=theta
    alphag(isave,:)=alpha
    psig(isave,:)=psi
    omegag(isave,:)=omega

    wbg(isave,:)=wb
    fxobsg(:,:,isave)=fxobs
    fxgridg(:,:,isave)=fxgrid
    yest=wb+sum(fxobs,2)+mu
    yestg(isave,:)=yest

    mug(isave,:)=mu
    mu0g(isave)=mu0
    kappag(isave)=kappa
    sigma2g(isave,:)=sigma2
    tmassg(isave)=tmass
    config(isave,:)=classind
    nclassg(isave)=nclass
    invlikeg(isave,:)=invlike

    if (verbose.eq.1) then
      if (mod(isave,ndisp).eq.0) then
        call cpu_time(itime)
        call sprint(isave,smcmc,itime-stime)
      end if
    end if
    isave=isave+1
  end if
end do
deallocate(phixobsfree,phixobsfreet,phixobs,phixgridfree,phixgrid)
pmetg=dble(pmet)/dble(smcmc*nskip)
if (verbose.eq.1) then
  do ifun=1,nfun
    if(fmodel(ifun).gt.1) then
      call aprint(ifun, pmetg(ifun))
    end if
  end do
end if

allocate(e_index(nobs))
mu0m=sum(mu0g)/dble(smcmc)
do imcmc=1,smcmc
  tmass=tmassg(imcmc)
  kappa=kappag(imcmc)
  mu0=mu0g(imcmc)
  nclass=nclassg(imcmc)
  allocate(muh(nclass),sigma2h(nclass),nclassh(nclass))
  classind=config(imcmc,:)
  mu=mug(imcmc,:)
  sigma2=sigma2g(imcmc,:)
  do iclass=1,nclass
    call which(classind.eq.iclass,nobs,e_index,nclassh(iclass))
    muh(iclass)=mu(e_index(1))
    sigma2h(iclass)=sigma2(e_index(1))
  end do
  muh=muh-mu0m
  mu0=mu0-mu0m

  call dp_update_predictive_dens(egrid,nclassh,muh,sigma2h,tmass,mu0,kappa,sigma2_r0,sigma2_s0,&
                                 nclass,nobs,ngrid,edens)
  edensg(imcmc,:)=edens
  deallocate(muh,sigma2h,nclassh)
end do
deallocate(e_index)
call rndend()

!=========================================================================================
contains

!======= Full conditionals ===============================================================
subroutine GetMCMC()
implicit none

!Internal arguments
real(8) :: gamrnd,rndunif,normrnd,ltnormrnd,rtnormrnd,tnormrnd,ltgamrnd,cdfnorm
real(8) :: resid(nobs),vni(nparw,nparw),vn(nparw,nparw),bn(nparw)
real(8) :: ck,rk(nobs),vni1(nbasis,nbasis),vn1(nbasis,nbasis),bn1(nbasis)
real(8) :: met_var_all_new(nfun),met_beta_new(nfun),testp
real(8) :: theta0_old,theta_old(nbasis),met_var0,met_var,met_std0
real(8) :: theta0_new,theta02_new,theta_new(nbasis),theta2_new(nbasis),thetanew(nbasis+1)
real(8) :: theta0_new_lnp,theta0_old_lnp,fxobs_new(nobs),fxgrid_new(nint+1)
real(8) :: met_stdS,met_varS,psi_old,psi_lnpold,psi_new,psi_lnpnew
real(8) :: omega_old,omega_lnpold,omega_new,omega_lnpnew
real(8) :: resid_old(nobs),resid_new(nobs),sse_old,sse_new,sold,snew
real(8) :: fx1(nobs),fxg1(nint+1),a_vni,a_vn,a_mn
real(8) :: th2gam(nbasis),sn,bup
real(8) :: ck1(nbasis),ckz(nbasis),u1(nbasis),u2,bnz(nbasis),bmin
real(8),allocatable :: bnz1(:)
integer :: z(nbasis)


if(ndimw.ge.1) then
  resid=yobs-sum(fxobs,2)-mu
  vni=matmul(wdatat,matmul(sigmai,wdata))
  vni=vni+beta_v0i
  call inverse(vni,nparw,vn)
  bn=matmul(vn,beta_v0im0+matmul(wdatat,resid/sigma2))

  call mvnrnd(bn,vn,nparw,beta)
  wb=matmul(wdata,beta)
else
  beta=0.d0
  wb=0.d0
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
resid=yobs-wb-sum(fxobs,2)-mu
ifree=1
irest=1
if (iflag_AM.eq.1) then
  icount_AM=icount_AM+1
end if
do ifun=1,nfun
  testp=0.d0
  rk=resid+fxobs(:,ifun)
  if (fmodel(ifun).eq.1) then
    vni1=matmul(phixobsfreet(:,:,ifree),matmul(sigmai,phixobsfree(:,:,ifree)))
    do k=1,nbasis
      vni1(k,k)=vni1(k,k)+1.d0/thv0(k,ifun)
    end do
    call inverse(vni1,nbasis,vn1)
    bn1=matmul(vn1,matmul(phixobsfreet(:,:,ifree),rk/sigma2))
    call mvnrnd(bn1,vn1,nbasis,theta(2:(nbasis+1),ifun))
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
    met_std0=dsqrt(met_var0)

    theta0_new=ltnormrnd(theta0_old,met_std0,0.d0)
    do k=1,nbasis
      theta_new(k)=normrnd(theta_old(k),dsqrt(met_var*thv0(k,ifun)))
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
    sse_new=sum(resid_new**2.d0/sigma2)

    resid_old=rk-fxobs(:,ifun)
    sse_old=sum(resid_old**2.d0/sigma2)

    snew=sum(theta2_new/thv0(:,ifun))
    sold=sum(theta2(:,ifun)/thv0(:,ifun))

    testp=testp- &
          (sse_new-sse_old)/(2.d0)- &
          (theta02_new)/(2.d0*theta0_v0)+ &
          (theta02(ifun))/(2.d0*theta0_v0)- &
          (snew-sold)/(2.d0)- &
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
      met_mAM(ifun)=met_mAM(ifun)+(met_var_all(ifun)-met_mAM(ifun))/dble(icount_AM)
      met_vAM(ifun)=(dble(icount_AM-1)/dble(icount_AM))*met_vAM(ifun) + &
                    ((met_var_all(ifun)-met_mAM(ifun))**2.d0)/dble(icount_AM)
    end if
    ck=metw*met_var_all(ifun)+(1.d0-metw)*met_mAM(ifun)
    met_beta_AM(ifun)=(met_alpha-1.d0)*ck

    irest=irest+1
  end if
  resid=yobs-wb-sum(fxobs,2)-mu
end do

do ifun=1,nfun
  if(fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
    resid=yobs-wb-sum(fxobs,2)-mu
    fx1=fxobs(:,ifun)-alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxg1=fxgrid(:,ifun)-alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))

    rk=resid+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))

    a_vni=sum((xobs(:,ifun)-xmid(ifun))**2.d0/sigma2)
    a_vni=a_vni+alpha_v0i
    a_vn=1.d0/a_vni
    a_mn=a_vn*(alpha_v0im0+dot_product((xobs(:,ifun)-xmid(ifun)),rk/sigma2))
    if (fpm(ifun).eq.1.d0) then
      alpha(ifun)=ltnormrnd(a_mn,dsqrt(a_vn),0.d0)
    else
      alpha(ifun)=rtnormrnd(a_mn,dsqrt(a_vn),0.d0)
    end if
    fxobs(:,ifun)=fx1+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxgrid(:,ifun)=fxg1+alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))
  else
    alpha(ifun)=0.d0
  end if
end do

do ifun=1,nfun
  th2gam=theta2(:,ifun)/gamvec(:,ifun)
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

  ck1=theta2(:,ifun)/(2.d0*tau2(ifun))

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


resid=yobs-wb-sum(fxobs,2)
call dp_update_configs_nogaps(resid,tmass,mu0,kappa,sigma2_r0,sigma2_s0,nobs,nclass,classind,mu,sigma2)

call dp_update_totalmass(tmass_a,tmass_b,nclass,nobs,tmass)

allocate(muh(nclass),sigma2h(nclass),nclassh(nclass))
call dp_update_classparams(resid,classind,mu0,kappa,sigma2_r0,sigma2_s0,nobs,nclass,muh,sigma2h,nclassh,mu,sigma2)
call diagvec(1.d0/sigma2,nobs,sigmai)

call dp_update_mu0(muh,sigma2h,kappa,nclass,mu0)

call dp_update_kappa(muh,sigma2h,mu0,kappa_r0,kappa_s0,nclass,kappa)

if(imcmc.gt.nblow .and. mod(imcmc,nskip).eq.0) then
  call dp_compute_invlike(resid,nclassh,muh,sigma2h,tmass,mu0,kappa,sigma2_r0,sigma2_s0,nclass,nobs,invlike)

end if
deallocate(muh,sigma2h,nclassh)

return
end subroutine GetMCMC


subroutine dp_update_configs_nogaps(resid,tmass,mu0,kappa,asig,bsig,nobs,nclass,classind,mu,sigma2)
implicit none

! Input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: resid(nobs),mu0,kappa,asig,bsig,tmass

! Output arguments
integer,intent(inout) :: classind(nobs),nclass
real(8), intent(inout) :: mu(nobs),sigma2(nobs)

! Internal arguments
integer :: iobs,iJ
integer :: ni,uns,si,J_i,s_i(nobs),sJ_i(nobs),re_s(nobs),s_n(nobs),s_index(nobs)
integer,allocatable :: nj_i(:)
real(8), allocatable :: logppf(:),loglik(:),logprob(:),prob(:),probn(:)
real(8), allocatable :: muStar(:),sigStar(:)
real(8) :: unspr
real(8) :: rndunif,normrnd,gamrnd,dnrm

s_n=classind
do iobs=1,nobs
  call rchkusr() ! check interrupt

  s_i=0
  s_i=s_n
  si=s_i(iobs)
  s_i(iobs)=0

  sJ_i=0
  call find_uniquei(s_i,nobs,sJ_i,J_i)
  ni=count(si.eq.s_i)

  uns=0
  unspr=(dble(J_i)-1.d0)/dble(J_i)
  if (ni.ne.0 .and. unspr.ge.rndunif()) uns=1
  if (uns.eq.1) then
    s_i(iobs)=si
    s_n=s_i
  else
    allocate(nj_i(J_i),muStar(J_i+1),sigStar(J_i+1))
    nj_i=0
    muStar=0.d0
    sigStar=0.d0
    re_s=0

    do iJ=1,J_i
      s_index=0
      call which(s_i.eq.sJ_i(iJ),nobs,s_index,nj_i(iJ))
      re_s(s_index(1:nj_i(iJ)))=iJ
      muStar(iJ)=mu(s_index(1))
      sigStar(iJ)=sigma2(s_index(1))
    end do

    sigStar(J_i+1)=1.d0/gamrnd(asig/2.d0,2.d0/bsig)
    muStar(J_i+1)=normrnd(mu0,dsqrt(kappa*sigStar(J_i+1)))

    allocate(loglik(J_i+1),logppf(J_i+1))
    allocate(logprob(J_i+1),prob(J_i+1),probn(J_i+1))
    loglik=0.d0
    logppf=0.d0
    logprob=0.d0
    prob=0.d0
    probn=0.d0

    logppf(1:J_i)=dlog(dble(nj_i))
    logppf(J_i+1)=dlog(tmass)-dlog(dble(J_i)+1.d0)

    do iJ=1,(J_i+1)
      loglik(iJ)=dnrm(resid(iobs),muStar(iJ),dsqrt(sigStar(iJ)),1)
    end do

    logprob=logppf+loglik
    prob=dexp(logprob-maxval(logprob))
    probn=prob/sum(prob)

    si=discrnd(J_i+1,probn)
    re_s(iobs)=si
    mu(iobs)=muStar(si)
    sigma2(iobs)=sigStar(si)
    s_n=re_s

    deallocate(nj_i,muStar,sigStar)
    deallocate(loglik,logppf,logprob,prob,probn)
  end if
end do
classind=s_n
nclass=maxval(classind)

return
end subroutine dp_update_configs_nogaps


subroutine dp_update_totalmass(parA,parB,nclass,nobs,totalmass)
implicit none

!input arguments
integer,intent(in) :: nclass,nobs
real(8), intent(in) :: parA,parB

!output argument
real(8),intent(inout) :: totalmass

!internal arguments
real(8) :: eta,wi,a,b
real(8) :: betarnd,gamrnd,rndunif

eta=betarnd(totalmass+1.d0,dble(nobs))

b=parB-dlog(eta)
wi=(parA+dble(nclass)-1.d0)/(parA+dble(nclass)-1.d0+dble(nobs)*(parB-dlog(eta)))
if(rndunif().lt.wi) then
  a=parA+dble(nclass)
else
  a=parA+dble(nclass)-1.d0
end if
totalmass=gamrnd(a,1.d0/b)

return
end subroutine dp_update_totalmass


subroutine dp_update_classparams(obs,ss,mu0,kappa,nu0,s0,nobs,nclass,class_m,class_v,nclassh,obs_m,obs_v)
implicit none

!input arguments
integer,intent(in) :: nobs,nclass,ss(nobs)
real(8), intent(in) :: obs(nobs),mu0,kappa,nu0,s0

!output arguments
real(8), intent(inout) :: obs_m(nobs),obs_v(nobs)
real(8), intent(out) :: class_m(nclass),class_v(nclass)
integer,intent(out) :: nclassh(nclass)

!internal arguments
integer :: iobs,h
real(8) :: mu_mh,kappah,var_rh,var_sh
real(8) :: sobsh,bobsh,sobshsq
real(8) :: gamrnd,normrnd

class_m=0.d0
class_v=0.d0
nclassh=0
do h=1,nclass
  call rchkusr()

  nclassh(h)=count(ss.eq.h)
  sobsh=0.d0
  sobshsq=0.d0
  do iobs=1,nobs
    if(ss(iobs).eq.h) then
      sobsh=sobsh+obs(iobs)
      sobshsq=sobshsq+obs(iobs)**2.d0
    end if
  end do
  bobsh=sobsh/dble(nclassh(h))
  kappah=1.d0/(dble(nclassh(h))+1.d0/kappa)
  var_rh=nu0+dble(nclassh(h))
  var_sh=s0+sobshsq-(sobsh**2.d0)/dble(nclassh(h))
  var_sh=var_sh+kappah*(dble(nclassh(h))/kappa)*(bobsh-mu0)**2.d0
  mu_mh=kappah*(mu0/kappa+sobsh)

  class_v(h)=1.d0/gamrnd(var_rh/2.d0,2.d0/var_sh)

  class_m(h)=normrnd(mu_mh,dsqrt(kappah*class_v(h)))

  do iobs=1,nobs
    if(ss(iobs).eq.h) then
      obs_m(iobs)=class_m(h)
      obs_v(iobs)=class_v(h)
    end if
  end do
end do

return
end subroutine dp_update_classparams


subroutine dp_update_mu0(class_m,class_v,kappa,nclass,mu0)
implicit none

!input arguments
integer,intent(in) :: nclass
real(8), intent(in) :: class_m(nclass),class_v(nclass),kappa

!output argument
real(8),intent(out) :: mu0

!internal arguments
real(8) :: mu0_m,mu0_v,mu0_s,normrnd

mu0_v=1.d0/(0.01d0+sum(1.d0/class_v)/kappa)
mu0_s=dsqrt(mu0_v)
mu0_m=mu0_v*sum(class_m/class_v)/kappa

mu0=normrnd(mu0_m,mu0_s)

return
end subroutine dp_update_mu0


subroutine dp_update_kappa(class_m,class_v,mu0,kappa_r0,kappa_s0,nclass,kappa)
implicit none

!input arguments
integer,intent(in) :: nclass
real(8), intent(in) :: class_m(nclass),class_v(nclass),mu0,kappa_r0,kappa_s0

!output argument
real(8),intent(out) :: kappa

!internal arguments
real(8) :: kappa_rn,kappa_sn,gamrnd

kappa_rn=kappa_r0+dble(nclass)
kappa_sn=kappa_s0+sum((class_m-mu0)**2.d0/class_v)

kappa=1.d0/gamrnd(kappa_rn/2.d0,2.d0/kappa_sn)

return
end subroutine dp_update_kappa


subroutine dp_compute_invlike(obs,nclassh,muh,sigma2h,tmass,mu0,kappa,nu0,s0,nclass,nobs,invlike)
implicit none

!input arguments
integer,intent(in) :: nclass,nobs,nclassh(nclass)
real(8), intent(in) :: obs(nobs),muh(nclass),sigma2h(nclass)
real(8), intent(in) :: tmass,mu0,kappa,nu0,s0

!output arguments
real(8),intent(out) :: invlike(nobs)

!internal arguments
integer :: iobs,iclass
real(8) :: mu_new,sigma2_new,ppf(nclass+1),like(nobs)
real(8) :: gamrnd,normrnd,dnrm

ppf=0.d0
ppf(1:nclass)=dble(nclassh)
ppf(nclass+1)=tmass
ppf=ppf/(tmass+dble(nobs))

like=0.d0
do iobs=1,nobs
  do iclass=1,nclass
    like(iobs)=like(iobs)+ppf(iclass)*dnrm(obs(iobs),muh(iclass),dsqrt(sigma2h(iclass)),0)
  end do
  sigma2_new=1.d0/gamrnd(nu0/2.d0,2.d0/s0)
  mu_new=normrnd(mu0,dsqrt(kappa*sigma2_new))
  like(iobs)=like(iobs)+ppf(nclass+1)*dnrm(obs(iobs),mu_new,dsqrt(sigma2_new),0)
end do

invlike=1.d0/like

return
end subroutine dp_compute_invlike


subroutine dp_update_predictive_dens(egrid,nclassh,muh,sigma2h,tmass,mu0,kappa,nu0,s0,nclass,&
                                     nobs,ngrid,edens)
implicit none

!input arguments
integer,intent(in) :: nclass,ngrid,nobs,nclassh(nclass)
real(8), intent(in) :: egrid(ngrid),muh(nclass),sigma2h(nclass)
real(8), intent(in) :: tmass,mu0,kappa,nu0,s0

!output arguments
real(8),intent(out) :: edens(ngrid)

!internal arguments
integer :: igrid,iclass
real(8) :: ppf(nclass+1),mu_class,mu_new,sigma2_new
real(8) :: normrnd,gamrnd,dnrm

ppf=0.d0
ppf(1:nclass)=dble(nclassh)
ppf(nclass+1)=tmass
ppf=ppf/(tmass+dble(nobs))

edens=0.d0
do igrid=1,ngrid
  do iclass=1,nclass
    mu_class=muh(iclass)
    edens(igrid)=edens(igrid)+ppf(iclass)*dnrm(egrid(igrid),mu_class,dsqrt(sigma2h(iclass)),0)
  end do
  sigma2_new=1.d0/gamrnd(nu0/2.d0,2.d0/s0)
  mu_new=normrnd(mu0,dsqrt(kappa*sigma2_new))
  edens(igrid)=edens(igrid)+ppf(nclass+1)*dnrm(egrid(igrid),mu_new,dsqrt(sigma2_new),0)
end do

return
end subroutine dp_update_predictive_dens

end subroutine bsaramdplocationscale
