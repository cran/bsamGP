subroutine bsaramgetlogg(fmodel,betag,sigma2g,thetag,tau2g,gamparg,alphag,psig,omegag,&
                         smcmc,nparw,nfun,nbasis,nExtremes,maxNext,iflagpsi,&
                         beta_mn,beta_covi,lndetbcov,sigma2_rn,sigma2_sn,&
                         theta_mn,theta_sn,theta0_lnpn,tau2_rn,tau2_sn,&
                         gamma_mn,gamma_vn,gamma_lnpn,alpha_mn,alpha_vn,alpha_lnpn,&
                         psi_mn,psi_vn,psi_lnpn,omega_mn,omega_vn,omega_lnpn,logg)
  use bsamTools
  implicit none

  !input arguments
  integer,intent(in) :: nparw,nfun,nbasis,nExtremes(nfun),maxNext,smcmc,fmodel(nfun),iflagpsi
  real(8), intent(in) :: betag(smcmc,nparw),beta_mn(nparw)
  real(8), intent(in) :: beta_covi(nparw,nparw),lndetbcov
  real(8), intent(in) :: sigma2g(smcmc),sigma2_rn,sigma2_sn
  real(8), intent(in) :: thetag(nbasis+1,nfun,smcmc),theta_mn(nbasis+1,nfun)
  real(8), intent(in) :: theta_sn(nbasis+1,nfun),theta0_lnpn(nfun)
  real(8), intent(in) :: tau2g(smcmc,nfun),tau2_rn(nfun),tau2_sn(nfun)
  real(8), intent(in) :: gamparg(smcmc,nfun),gamma_mn(nfun),gamma_vn(nfun),gamma_lnpn(nfun)
  real(8), intent(in) :: alphag(smcmc,nfun),alpha_mn(nfun),alpha_vn(nfun),alpha_lnpn(nfun)
  real(8), intent(in) :: psig(smcmc,maxNext,nfun),psi_mn(maxNext,nfun),psi_vn(maxNext,nfun),psi_lnpn(nfun)
  real(8), intent(in) :: omegag(smcmc,maxNext,nfun),omega_mn(maxNext,nfun),omega_vn(maxNext,nfun),omega_lnpn(nfun)

  !output argument
  real(8),intent(out) :: logg(smcmc)

  !internal arguments
  integer :: imcmc,ifun,k
  real(8) :: beta(nparw),sigma2,theta(nbasis+1,nfun),tau2(nfun),gampar(nfun),alpha(nfun)
  real(8) :: psi(maxNext,nfun),omega(maxNext,nfun),residbeta(nparw),residtheta(nbasis+1),logIlik

  do imcmc=1,smcmc
    beta=betag(imcmc,:)
    sigma2=sigma2g(imcmc)
    theta=thetag(:,:,imcmc)
    tau2=tau2g(imcmc,:)
    gampar=gamparg(imcmc,:)
    alpha=alphag(imcmc,:)
    psi=psig(imcmc,:,:)
    omega=omegag(imcmc,:,:)

    logIlik=0.d0

    residbeta=beta-beta_mn
    logIlik=logIlik-dot_product(residbeta,matmul(beta_covi,residbeta))/2.d0- &
            dble(nparw)*dlog(2.d0*PI)/2.d0-lndetbcov/2.d0

    logIlik=logIlik + LogfIG(sigma2,sigma2_rn,sigma2_sn)

    do ifun=1,nfun
      if(fmodel(ifun).eq.1) then
        residtheta(2:(nbasis+1))=(theta(2:(nbasis+1),ifun)- &
                                  theta_mn(2:(nbasis+1),ifun))/theta_sn(2:(nbasis+1),ifun)
        logIlik=logIlik-sum(residtheta(2:(nbasis+1))**2.d0)/2.d0- &
                dlog(2.d0*PI)/2.d0-sum(dlog(theta_sn(2:(nbasis+1),ifun)))
      else
        residtheta=(theta(:,ifun)-theta_mn(:,ifun))/theta_sn(:,ifun)
        logIlik=logIlik-sum(residtheta**2.d0)/2.d0- &
                dble(nbasis+1)*dlog(2.d0*PI)/2.d0-sum(dlog(theta_sn(:,ifun)))

        logIlik=logIlik-theta0_lnpn(ifun)
      end if

      logIlik=logIlik+LogfIG(tau2(ifun),tau2_rn(ifun),tau2_sn(ifun))

      logIlik=logIlik-((gampar(ifun)-gamma_mn(ifun))**2.d0)/(2.d0*gamma_vn(ifun))- &
              dlog(2.d0*PI*gamma_vn(ifun))/2.d0-gamma_lnpn(ifun)

      if (fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
        logIlik=logIlik-((alpha(ifun)-alpha_mn(ifun))**2.d0)/(2.d0*alpha_vn(ifun))- &
                dlog(2.d0*PI*alpha_vn(ifun))-alpha_lnpn(ifun)
      end if

      if (fmodel(ifun).ge.5 .and. fmodel(ifun).le.7) then
        if(iflagpsi.eq.1) then
          logIlik=logIlik-((psi(1,ifun)-psi_mn(1,ifun))**2.d0)/(2.d0*psi_vn(1,ifun))- &
                  dlog(2.d0*PI*psi_vn(1,ifun))/2.d0-psi_lnpn(ifun)
        end if
        logIlik=logIlik-((omega(1,ifun)-omega_mn(1,ifun))**2.d0)/(2.d0*omega_vn(1,ifun))- &
                dlog(2.d0*PI*omega_vn(1,ifun))/2.d0-omega_lnpn(ifun)
      end if
      if (fmodel(ifun).eq.8) then
        do k=1,nExtremes(ifun)
          if(iflagpsi.eq.1) then
            logIlik=logIlik-((psi(k,ifun)-psi_mn(k,ifun))**2.d0)/(2.d0*psi_vn(k,ifun))- &
                    dlog(2.d0*PI*psi_vn(k,ifun))/2.d0
          end if
          logIlik=logIlik-((omega(k,ifun)-omega_mn(k,ifun))**2.d0)/(2.d0*omega_vn(k,ifun))- &
                  dlog(2.d0*PI*omega_vn(k,ifun))/2.d0
        end do
        logIlik=logIlik-psi_lnpn(ifun)-omega_lnpn(ifun)
      end if
    end do
    logg(imcmc)=logIlik
  end do

  return
end subroutine bsaramgetlogg
