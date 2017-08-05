subroutine gbsaramgetlogg(fmodel,betag,thetag,tau2g,gamparg,alphag,psig,omegag,&
                          smcmc,nparw,nfun,nbasis,iflagpsi,beta_mn,beta_covi,lndetbcov,&
                          theta_mn,theta_sn,theta0_lnpn,tau2_rn,tau2_sn,&
                          gamma_mn,gamma_vn,gamma_lnpn,alpha_mn,alpha_vn,alpha_lnpn,&
                          psi_mn,psi_vn,psi_lnpn,omega_mn,omega_vn,omega_lnpn,logg)
use bsamTools
implicit none

!input arguments
integer,intent(in) :: nparw,nfun,nbasis,smcmc,fmodel(nfun),iflagpsi
real(8), intent(in) :: betag(smcmc,nparw),beta_mn(nparw)
real(8), intent(in) :: beta_covi(nparw,nparw),lndetbcov
real(8), intent(in) :: thetag(nbasis+1,nfun,smcmc),theta_mn(nbasis+1,nfun)
real(8), intent(in) :: theta_sn(nbasis+1,nfun),theta0_lnpn(nfun)
real(8), intent(in) :: tau2g(smcmc,nfun),tau2_rn(nfun),tau2_sn(nfun)
real(8), intent(in) :: gamparg(smcmc,nfun),gamma_mn(nfun),gamma_vn(nfun),gamma_lnpn(nfun)
real(8), intent(in) :: alphag(smcmc,nfun),alpha_mn(nfun),alpha_vn(nfun),alpha_lnpn(nfun)
real(8), intent(in) :: psig(smcmc,nfun),psi_mn(nfun),psi_vn(nfun),psi_lnpn(nfun)
real(8), intent(in) :: omegag(smcmc,nfun),omega_mn(nfun),omega_vn(nfun),omega_lnpn(nfun)

!output argument
real(8),intent(out) :: logg(smcmc)

!internal arguments
integer :: imcmc,ifun
real(8) :: beta(nparw),theta(nbasis+1,nfun),tau2(nfun),gampar(nfun),alpha(nfun)
real(8) :: psi(nfun),omega(nfun),residbeta(nparw),residtheta(nbasis+1),logIlik

do imcmc=1,smcmc
  call rchkusr() ! check interrupt

  beta=betag(imcmc,:)
  theta=thetag(:,:,imcmc)
  tau2=tau2g(imcmc,:)
  gampar=gamparg(imcmc,:)
  alpha=alphag(imcmc,:)
  psi=psig(imcmc,:)
  omega=omegag(imcmc,:)

  logIlik=0.d0

  residbeta=beta-beta_mn
  logIlik=logIlik-dot_product(residbeta,matmul(beta_covi,residbeta))/2.d0- &
          dble(nparw)*dlog(2.d0*PI)/2.d0-lndetbcov/2.d0

  do ifun=1,nfun
    residtheta=(theta(:,ifun)-theta_mn(:,ifun))/theta_sn(:,ifun)
    if(fmodel(ifun).eq.1) then
      logIlik=logIlik-sum(residtheta(2:(nbasis+1))**2.d0)/2.d0- &
              dlog(2.d0*PI)/2.d0-sum(dlog(theta_sn(2:(nbasis+1),ifun)))
    else
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

    if (fmodel(ifun).ge.5) then
      if(iflagpsi.eq.1) then
        logIlik=logIlik-((psi(ifun)-psi_mn(ifun))**2.d0)/(2.d0*psi_vn(ifun))- &
                 dlog(2.d0*PI*psi_vn(ifun))/2.d0-psi_lnpn(ifun)
      end if
      logIlik=logIlik-((omega(ifun)-omega_mn(ifun))**2.d0)/(2.d0*omega_vn(ifun))- &
              dlog(2.d0*PI*omega_vn(ifun))/2.d0-omega_lnpn(ifun)
    end if
  end do
  logg(imcmc)=logIlik
end do

return
end subroutine gbsaramgetlogg
