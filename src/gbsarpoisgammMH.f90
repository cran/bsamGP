subroutine gbsarpoisgammMH(verbose,yint,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,&
                           theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,&
                           alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,&
                           kappa_m0,kappa_v0,iflagprior,iflagpsi,&
                           maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,&
                           zetag,tau2g,gammag,thetag,betag,alphag,psig,omegag,kappag,&
                           fxgridg,fxobsg,muhatg,wbg,loglikeg,logpriorg,imodmetg,pmetg)
use ToolsRfunf
use gbsamTools
implicit none

!input arguments
integer,intent(in) :: nobs,nparw,nfun,nbasis,nint,fmodel(nfun)
integer,intent(in) :: maxmodmet,nblow0,nblow,smcmc,nskip,ndisp
integer,intent(in) :: yint(nobs),iflagprior,iflagpsi,verbose
real(8), intent(in) :: wdata(nobs,nparw),xobs(nobs,nfun),fpm(nfun)
real(8), intent(in) :: theta0_m0,theta0_s0
real(8), intent(in) :: tau2_m0,tau2_v0,w0,beta_m0(nparw),beta_v0(nparw,nparw)
real(8), intent(in) :: alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed
real(8), intent(in) :: omega_m0(nfun),omega_s0(nfun),kappa_m0,kappa_v0

!output arguments
integer,intent(out) :: imodmetg
real(8),intent(out) :: zetag(smcmc,nfun),tau2g(smcmc,nfun),gammag(smcmc,nfun)
real(8),intent(out) :: thetag(nbasis+1,nfun,smcmc),betag(smcmc,nparw)
real(8),intent(out) :: alphag(smcmc,nfun),psig(smcmc,nfun),omegag(smcmc,nfun)
real(8),intent(out) :: fxobsg(nobs,nfun,smcmc),fxgridg(nint+1,nfun,smcmc)
real(8),intent(out) :: wbg(smcmc,nobs),muhatg(smcmc,nobs),kappag(smcmc)
real(8),intent(out) :: loglikeg(smcmc),logpriorg(smcmc),pmetg(nfun)

!internal arguments
real(8)  :: stime,itime
integer :: imcmc,isave,nmcmc

integer :: iloop,ifun,iobs,intsimpfacts(nint+1),xinxgrid(nobs,nfun)
real(8) :: xobsint(nobs,nfun),xidelta(nobs,nfun)

real(8) :: yobs(nobs),xmin(nfun),xmax(nfun),xrange(nfun),xmid(nfun),xdelta(nfun)
real(8) :: wtw(nparw,nparw),wtwi(nparw,nparw),wdatat(nparw,nobs)
real(8) :: xgrid(nint+1,nfun),xbar(nfun),xtx(nfun)

integer :: k,nfunconstraint
real(8) :: kall(nbasis),kall0(nbasis+1),kbar

real(8) :: phixobs(nobs,nbasis+1,nfun),phixobst(nbasis+1,nobs,nfun)
real(8) :: phixgrid(nint+1,nbasis+1,nfun),phi2(nbasis,nbasis,nfun)

real(8) :: theta0_v0,theta0_v0i,theta0_v0im0,theta0_lnv0
real(8) :: tau2_r0,tau2_s0,tau2_rn,tau2_u0,wk
real(8) :: beta_v0i(nparw,nparw),beta_v0im0(nparw),beta_lnv0
real(8) :: alpha_v0,alpha_v0i,alpha_v0im0,alpha_lnv0
real(8) :: psi_v0i,psi_v0,psi_lnp0,psi_lnv0
real(8) :: omega_v0i(nfun),omega_v0(nfun),omega_lnp0(nfun),omega_lnv0(nfun)
real(8) :: kappa_r0,kappa_s0

integer :: imodmet,pmet(nfun),iflag_AM,icount_AM,pok,ibeta_AM
integer :: nbatch,ibatch,icount_alpha,a_nacc(nfun) ! alpha
integer :: kappa_nacc
real(8) :: metw,met_alpha,met_omega_alpha,met_psi_alpha
real(8) :: metm(nfun),mets(nfun),metv(nfun),met_beta(nfun),met_beta_AM(nfun)
real(8) :: met_mAM(nfun),met_vAM(nfun),met_omega_m(nfun),met_omega_beta(nfun)
real(8) :: met_psi_m(nfun),met_psi_beta(nfun),met_var_all(nfun),ck
real(8) :: Iw(nparw,nparw),beta_mAM(nparw,1),beta_m1AM(nparw,1),beta_vAM(nparw,nparw)
real(8) :: a_lsd(nfun),kappa_lsd
real(8),allocatable :: betas(:,:)

real(8) :: kappa,tau(nfun),tau2(nfun),tau2i(nfun),gampar(nfun),lngampar(nfun)
real(8) :: theta(nbasis+1,nfun),theta02(nfun),theta2(nbasis,nfun),beta(nparw),wb(nobs)
real(8) :: theta0(nfun),gamvec0(nbasis+1,nfun),thv0(nbasis,nfun),thv00(nbasis+1,nfun)
real(8) :: psi(nfun),omega(nfun),zeta(nfun),alpha(nfun),gamvec(nbasis,nfun)
real(8) :: fxobs(nobs,nfun),fxgrid(nint+1,nfun),yest(nobs),lambda(nobs)

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

call diag(1.d0,nparw,Iw)

nbatch=100
ibatch=0
icount_alpha=0
a_nacc=0
a_lsd=0.d0

kappa_nacc=0
kappa_lsd=0.d0

yobs=dble(yint)
wdatat=transpose(wdata)
wtw=matmul(wdatat,wdata)
call inverse(wtw,nparw,wtwi)
xbar=sum(xobs,1)/dble(nobs)
do ifun=1,nfun
  xmin(ifun)=minval(xobs(:,ifun))
  xmax(ifun)=maxval(xobs(:,ifun))
  xmid(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
  xrange(ifun)=xmax(ifun)-xmin(ifun)
  xdelta(ifun)=(xrange(ifun))/dble(nint)
  xtx(ifun)=sum((xobs(:,ifun)-xmid(ifun))**2.d0)
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

kall=(/ (dble(k),k=1,nbasis) /)
kall0=(/ (dble(k),k=0,nbasis) /)
kbar=sum(kall0)/dble(nbasis+1)

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

kappa_r0=2.d0*(2.d0+kappa_m0**2.d0/kappa_v0)
kappa_s0=kappa_m0*(kappa_r0-2.d0)

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
  if(fmodel(ifun).eq.1) then
    theta(1,ifun)=0.d0
  else
    theta(1,ifun)=0.1d0
  end if
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
kappa=kappa_m0

nfunconstraint=count(fmodel.gt.1)
do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobs(:,2:(nbasis+1),ifun))
    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgrid(:,2:(nbasis+1),ifun))
    phixobs(:,1,ifun)=0.d0
    phixgrid(:,1,ifun)=0.d0
    phixobst(:,:,ifun)=transpose(phixobs(:,:,ifun))
    phi2(:,:,ifun)=matmul(phixobst(2:(nbasis+1),:,ifun),phixobs(:,2:(nbasis+1),ifun))

    call GetFreef(theta(2:(nbasis+1),ifun),phixobs(:,2:(nbasis+1),ifun),&
                  phixgrid(:,2:(nbasis+1),ifun),nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))

  else
    phixobst(:,:,ifun)=0.d0
    phi2(:,:,ifun)=0.d0

    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis,phixobs(:,2:(nbasis+1),ifun))
    call ConstFun(xobs(:,ifun),xrange(ifun),nobs,phixobs(:,1,ifun))

    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis,phixgrid(:,2:(nbasis+1),ifun))
    call ConstFun(xgrid(:,ifun),xrange(ifun),nint+1,phixgrid(:,1,ifun))

    if (fmodel(ifun).eq.2) then
      call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun), &
                  xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun), &
                  intsimpfacts,nbasis+1,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
    else if (fmodel(ifun).eq.3) then
      call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
                      xmid(ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),xdelta(ifun),&
                      xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
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

if(maxval(fmodel).gt.1) then
  ibeta_AM=0
  if (verbose.eq.1) then
    call piprint()
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
    end if
  end do
end if

ibeta_AM=1
isave=1
nmcmc=nblow+nskip*smcmc
allocate(betas(nmcmc+1,nparw))
betas=0.d0
beta_mAM(:,1)=beta
do imcmc=1,nmcmc
  if(imcmc.eq.1) then
    if (verbose.eq.1) then
      call biprint()
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
      call miprint()
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
    kappag(isave)=kappa

    wbg(isave,:)=wb
    fxobsg(:,:,isave)=fxobs
    fxgridg(:,:,isave)=fxgrid
    yest=wb+sum(fxobs,2)
    muhatg(isave,:)=dexp(yest)

    loglikeg(isave)=GetLogLik_NegBin(yobs,dexp(yest),kappa,nobs)

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
pmetg=dble(pmet)/dble(smcmc*nskip)
if (verbose.eq.1) then
  do ifun=1,nfun
    if(fmodel(ifun).gt.1) then
      call aprint(ifun, pmetg(ifun))
    end if
  end do
end if
deallocate(betas)
call rndend()

!=========================================================================================
contains

!======= Full conditionals ===============================================================
subroutine GetMCMC()
implicit none

!Internal arguments
real(8) :: gamrnd,rndunif,normrnd,ltnormrnd,tnormrnd,ltgamrnd
real(8) :: lambda_r0(nobs),lambda_s0(nobs)
real(8) :: beta_new(nparw),beta_old(nparw),loglik_new,loglik_old
real(8) :: wb_old(nobs),wb_new(nobs),beta_vn(nparw,nparw),betat(1,nparw),beta_testp
real(8) :: rk(nobs)
real(8) :: met_var_all_new(nfun),met_beta_new(nfun),testp
real(8) :: theta0_old,theta_old(nbasis),met_var0,met_var,met_std0
real(8) :: theta0_new,theta02_new,theta_new(nbasis),theta2_new(nbasis),thetanew(nbasis+1)
real(8) :: theta0_new_lnp,theta0_old_lnp,fxobs_new(nobs),fxgrid_new(nint+1)
real(8) :: met_stdS,met_varS,psi_old,psi_lnpold,psi_new,psi_lnpnew
real(8) :: omega_old,omega_lnpold,omega_new,omega_lnpnew,sold,snew
real(8) :: fx1(nobs),fxg1(nint+1),alpha_old,alpha_new,alpha_testp
real(8) :: kappa_old,kappa_new,kappa_testp
real(8) :: th2gam(nbasis),sn,bup
real(8) :: ck1(nbasis),ckz(nbasis),u1(nbasis),u2,bnz(nbasis),bmin
real(8),allocatable :: bnz1(:)
integer :: z(nbasis)

lambda_r0=yobs+kappa
lambda_s0=kappa/dexp(wb+sum(fxobs,2))+1.d0
do iobs=1,nobs
  lambda(iobs)=gamrnd(lambda_r0(iobs),1.d0/lambda_s0(iobs))
end do

beta_old=beta
wb_old=wb
loglik_old=GetLogLik_Gamma(lambda,wb_old+sum(fxobs,2),kappa,nobs)

if(ibeta_AM.eq.0) then
  beta_vn=(0.01d0)*Iw/dble(nparw)
else if(ibeta_AM.eq.1 .and. imcmc.le.nparw*2) then
  beta_vn=(0.01d0)*Iw/dble(nparw)
else if(ibeta_AM.eq.1 .and. imcmc.gt.nparw*2) then
  if(rndunif().le.0.05d0) then
    beta_vn=(0.01d0)*Iw/dble(nparw)
  else
    beta_vn=((2.38d0)**2.d0)*beta_vAM/dble(nparw)
  end if
end if
call mvnrnd(beta_old,beta_vn,nparw,beta_new)
wb_new=matmul(wdata,beta_new)
loglik_new=GetLogLik_Gamma(lambda,wb_new+sum(fxobs,2),kappa,nobs)

beta_testp=loglik_new+mvnpdf(beta_new,beta_m0,beta_v0,nparw,.true.)- &
loglik_old-mvnpdf(beta_old,beta_m0,beta_v0,nparw,.true.)
if(dlog(rndunif()).le.beta_testp) then
  beta=beta_new
  wb=matmul(wdata,beta)
end if
if(ibeta_AM.eq.1) betas(imcmc+1,:)=beta

if(ibeta_AM.eq.1) then
  beta_m1AM=beta_mAM
  beta_mAM(:,1)=(dble(imcmc)*beta_m1AM(:,1)+beta)/dble(imcmc+1)
  if(imcmc.eq.nparw*2) call covariance(betas(1:(imcmc+1),:),imcmc+1,nparw,beta_vAM)
  if(imcmc.gt.nparw*2) then
    betat(1,:)=beta
    beta_vAM=(dble(imcmc-1)/dble(imcmc))*beta_vAM+ &
             (dble(imcmc)*matmul(beta_m1AM,transpose(beta_m1AM))- &
             dble(imcmc-1)*matmul(beta_mAM,transpose(beta_mAM))+ &
             matmul(transpose(betat),betat))/dble(imcmc)
  end if
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
if (iflag_AM.eq.1) then
  icount_AM=icount_AM+1
end if
do ifun=1,nfun
  testp=0.d0
  rk=wb+sum(fxobs,2)-fxobs(:,ifun)

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
  if (fmodel(ifun).eq.1) then
    call GetFreef(thetanew(2:(nbasis+1)),phixobs(:,2:(nbasis+1),ifun),phixgrid(:,2:(nbasis+1),ifun),&
                  nbasis,nobs,nint+1,fxobs_new,fxgrid_new)
  else if (fmodel(ifun).eq.2) then
    call GetUpf(fpm(ifun),thetanew,phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                intsimpfacts,nbasis+1,nobs,nint+1,fxobs_new,fxgrid_new)
  else if (fmodel(ifun).eq.3) then
    call GetConvexf(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
                    xmid(ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                    xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                    intsimpfacts,nbasis+1,nobs,nint+1,fxobs_new,fxgrid_new)
  else if (fmodel(ifun).eq.4) then
    call GetConcavef(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
                     xmid(ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                     xdelta(ifun),xrange(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                     intsimpfacts,nbasis+1,nobs,nint+1,fxobs_new,fxgrid_new)
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
                 xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun), &
                 xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),&
                 xmid(ifun),intsimpfacts,nbasis+1,nobs,nint+1,fxobs_new,fxgrid_new)
    else if (fmodel(ifun).eq.6) then
      call GetRotateSf(fpm(ifun),omega_new,psi_new,alpha(ifun),thetanew, &
                       xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                       xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),&
                       xrange(ifun),xmid(ifun),intsimpfacts,nbasis+1,nobs,nint+1,&
                       fxobs_new,fxgrid_new)
    else
      call GetUf(fpm(ifun),omega_new,psi_new,thetanew,xobs(:,ifun), &
                 xgrid(:,ifun),phixobs(:,:,ifun),phixgrid(:,:,ifun),&
                 xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),&
                 intsimpfacts,nbasis+1,nobs,nint+1,fxobs_new,fxgrid_new)
    end if
  end if

  loglik_new=GetLogLik_Gamma(lambda,rk+fxobs_new,kappa,nobs)
  loglik_old=GetLogLik_Gamma(lambda,rk+fxobs(:,ifun),kappa,nobs)

  snew=sum(theta2_new/thv0(:,ifun))
  sold=sum(theta2(:,ifun)/thv0(:,ifun))

  if(fmodel(ifun).eq.1) then
    testp=testp+ &
          (loglik_new-loglik_old)- &
          (snew-sold)/(2.d0)- &
          ((met_alpha+1.d0)*dlog(met_var_all(ifun)))-&
          (met_beta_new(ifun)/met_var_all(ifun))+ &
          ((met_alpha+1.d0)*dlog(met_var_all_new(ifun)))+&
          (met_beta_AM(ifun)/met_var_all_new(ifun))
  else
    testp=testp+ &
          (loglik_new-loglik_old)- &
          (theta02_new)/(2.d0*theta0_v0)+ &
          (theta02(ifun))/(2.d0*theta0_v0)- &
          (snew-sold)/(2.d0)- &
          theta0_old_lnp+theta0_new_lnp- &
          ((met_alpha+1.d0)*dlog(met_var_all(ifun)))-&
          (met_beta_new(ifun)/met_var_all(ifun))+ &
          ((met_alpha+1.d0)*dlog(met_var_all_new(ifun)))+&
          (met_beta_AM(ifun)/met_var_all_new(ifun))
  end if

  if(dlog(rndunif()).le.testp) then
    if(fmodel(ifun).eq.1) then
      theta0(ifun)=0.d0
      theta02(ifun)=0.d0
      theta(1,ifun)=0.d0
    else
      theta0(ifun)=theta0_new
      theta02(ifun)=theta02_new
      theta(1,ifun)=theta0_new
    end if
    theta(2:(nbasis+1),ifun)=theta_new
    theta2(:,ifun)=theta2_new
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
end do

icount_alpha=icount_alpha+1
if(mod(icount_alpha,nbatch).eq.0) then
  ibatch=ibatch+1
  do ifun=1,nfun
    if(fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
      if(dble(a_nacc(ifun))/dble(icount_alpha).lt.0.44d0) then
        a_lsd(ifun)=a_lsd(ifun)-min(0.01d0,1.d0/dsqrt(dble(ibatch)))
      else if(dble(a_nacc(ifun))/dble(icount_alpha).gt.0.44d0) then
        a_lsd(ifun)=a_lsd(ifun)+min(0.01d0,1.d0/dsqrt(dble(ibatch)))
      end if
    end if
  end do
end if

do ifun=1,nfun
  rk=wb+sum(fxobs,2)
  if(fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
    fx1=fxobs(:,ifun)-alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxg1=fxgrid(:,ifun)-alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))

    rk=rk-alpha(ifun)*(xobs(:,ifun)-xmid(ifun))

    alpha_old=alpha(ifun)
    loglik_old=GetLogLik_Gamma(lambda,rk+alpha_old*(xobs(:,ifun)-xmid(ifun)),kappa,nobs)
    sold=(alpha_old-alpha_m0)**2.d0

    alpha_new=normrnd(alpha_old,dexp(a_lsd(ifun)))
    loglik_new=GetLogLik_Gamma(lambda,rk+alpha_new*(xobs(:,ifun)-xmid(ifun)),kappa,nobs)
    snew=(alpha_new-alpha_m0)**2.d0

    alpha_testp=loglik_new-loglik_old- &
                (snew-sold)/(2.d0*alpha_v0)
    if(dlog(rndunif()).le.alpha_testp) then
      if(fpm(ifun).gt.0.d0 .and. alpha_new.gt.0.d0) then
        alpha(ifun)=alpha_new
        a_nacc(ifun)=a_nacc(ifun)+1
      else if(fpm(ifun).lt.0.d0 .and. alpha_new.lt.0.d0) then
        alpha(ifun)=alpha_new
        a_nacc(ifun)=a_nacc(ifun)+1
      end if
    end if

    fxobs(:,ifun)=fx1+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
    fxgrid(:,ifun)=fxg1+alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))
  else
    alpha(ifun)=0.d0
  end if
end do

if(mod(icount_alpha,nbatch).eq.0) then
  if(dble(kappa_nacc)/dble(icount_alpha).lt.0.44d0) then
    kappa_lsd=kappa_lsd-min(0.01d0,1.d0/dsqrt(dble(ibatch)))
  else if(dble(kappa_nacc)/dble(icount_alpha).gt.0.44d0) then
    kappa_lsd=kappa_lsd+min(0.01d0,1.d0/dsqrt(dble(ibatch)))
  end if
end if

kappa_old=kappa
kappa_new=normrnd(kappa_old,dexp(kappa_lsd))

kappa_testp=GetLogPost_kappa(kappa_new,lambda,wb+sum(fxobs,2),nobs,kappa_r0,kappa_s0)- &
            GetLogPost_kappa(kappa_old,lambda,wb+sum(fxobs,2),nobs,kappa_r0,kappa_s0)

if(dlog(rndunif()).le.kappa_testp) then
  kappa=kappa_new
  kappa_nacc=kappa_nacc+1
end if

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

return
end subroutine GetMCMC

real(8) function GetLogLik_Gamma(lambda,logmu,kappa,nobs)
implicit none

!input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: lambda(nobs),logmu(nobs),kappa

GetLogLik_Gamma=-kappa*(sum(logmu)+sum(lambda*dexp(-logmu)))

return
end function GetLogLik_Gamma

real(8) function GetLogPost_kappa(kappa,lambda,Xb,nobs,r0,s0)
implicit none

! input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: kappa,lambda(nobs),Xb(nobs),r0,s0

! internal arguments
real(8) :: gammaln

GetLogPost_kappa=(dble(nobs)*kappa+r0-1.d0)*dlog(kappa)-dble(nobs)*gammaln(kappa)-&
                 kappa*(sum(Xb))+(kappa-1.d0)*sum(dlog(lambda))-kappa*(s0+sum(lambda*dexp(-Xb)))

return
end function GetLogPost_kappa

real(8) function GetLogLik_NegBin(y,mu,kappa,nobs)
implicit none

!input arguments
integer,intent(in) :: nobs
real(8), intent(in) :: y(nobs),mu(nobs),kappa

!internal argument
integer :: iobs
real(8) :: const(nobs),loglik,gammaln

do iobs=1,nobs
  const(iobs)=gammaln(y(iobs)+kappa)-gammaln(y(iobs)+1.d0)
end do
loglik=-sum(y*dlog(kappa/mu+1.d0))-kappa*sum(dlog(1.d0+mu/kappa))- &
       dble(nobs)*gammaln(kappa)+sum(const)
GetLogLik_NegBin=loglik

return
end function GetLogLik_NegBin


real(8) function GetLogPrior()
implicit none

!Internal argument
real(8) :: residb(nparw),thvar(nbasis+1),thetasq(nbasis+1,nfun)
real(8) :: sse,theta0_lnp0,alpha_lnp0,lnpriorf
real(8) :: cdfnorm,dgamm

lnpriorf=0.d0
thetasq=theta**2.d0

residb=beta-beta_m0
lnpriorf=lnpriorf-dot_product(residb,matmul(beta_v0i,residb))/(2.d0)- &
         dble(nparw)*dlog(2.d0*PI)/2.d0-beta_lnv0/2.d0

lnpriorf=lnpriorf+dgamm(kappa,kappa_r0,1.d0/kappa_s0,1)

do ifun=1,nfun
  if (fmodel(ifun).eq.1) then
    thvar(2:(nbasis+1))=tau2(ifun)*gamvec(:,ifun)
    sse=sum(thetasq(2:(nbasis+1),ifun)/thvar(2:(nbasis+1)))
    lnpriorf=lnpriorf-sse/2.d0-dble(nbasis)*dlog(2.d0*PI)/2.d0- &
             sum(dlog(thvar(2:(nbasis+1))))/2.d0
  else
    thvar(1)=theta0_v0
    thvar(2:(nbasis+1))=tau2(ifun)*gamvec(:,ifun)
    sse=sum(thetasq(:,ifun)/thvar)
    theta0_lnp0=cdfnorm(-theta0_m0/(theta0_s0),0.d0,1.d0,0,1)
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
    alpha_lnp0=cdfnorm(-alpha_m0/alpha_s0,0.d0,1.d0,0,1)
    lnpriorf=lnpriorf-((alpha(ifun)-alpha_m0)**2.d0)/(2.d0*alpha_v0)- &
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


end subroutine gbsarpoisgammMH
