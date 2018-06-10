#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/******************************************************************************************************
 Random number generators
 ******************************************************************************************************/

void F77_SUB(rndstart)(void)
{
    GetRNGstate();
}

void F77_SUB(rndend)(void)
{
    PutRNGstate();
}

double F77_SUB(rndunif)(void)
{
    return unif_rand();
}

double F77_SUB(rndnorm)(void)
{
    return norm_rand();
}

double F77_SUB(normrnd)(double *mu, double *sigma)
{
    return rnorm(*mu, *sigma);
}

double F77_SUB(rtnormrnd)(double *mu, double *sigma, double *up)
{
    double u,pcut,x;

    if ((*sigma)==0.0){ /* degenerate sigma=0 */
        return ((*mu)<(*up)) ? (*mu) : *(up);
    }
    pcut=pnorm(*up,*mu,*sigma,1,0);
    if (pcut < 0.0001)
        return (*up)-0.0001*(*sigma);
    u=unif_rand();
    u=u*pcut;
    x=qnorm(u,*mu,*sigma,1,0);

    return x;
}

double F77_SUB(ltnormrnd)(double *mu, double *sigma, double *low)
{
    double u,pcut,x;

    if ((*sigma)==0.0){ /* degenerate sigma=0 */
        return ((*mu)>(*low)) ? (*mu) : (*low);
    }

    pcut=pnorm(*low,*mu,*sigma,1,0);
    if (pcut>0.9999)
        return (*low)+0.0001*(*sigma);
    u=unif_rand();
    u=pcut+(1.0-pcut)*u;
    x=qnorm(u,*mu,*sigma,1,0);

    return x;
}

double F77_SUB(tnormrnd)(double *mu, double *sigma, double *low, double *up)
{
    double u, pleft, pright, y;

    pleft=pnorm(*low,*mu,*sigma,1,0);
    pright=pnorm(*up,*mu,*sigma,1,0);
    if (pleft>0.9999)
        return (*low)+0.0001*fmax2((*up)-(*low),*sigma);
    if (pright<0.0001)
        return (*up)-0.0001*fmax2((*up)-(*low),*sigma);
    u=unif_rand();
    u=pleft+(pright-pleft)*u;
    y=qnorm(u,*mu,*sigma,1,0);

    return y;
}

double F77_SUB(gamrnd)(double *shape, double *scale)
{
    return rgamma(*shape, *scale);
}

double F77_SUB(rtgamrnd)(double *shape, double *scale, double *up)
{
    double gup,p,x;
    gup=pgamma(*up,*shape,*scale,1,0);
    if(gup>=1.0) {
        gup=0.99999;
    }
    if(gup<=0.0) {
        gup=0.00001;
    }
    p=unif_rand()*gup;
    x=qgamma(p,*shape,*scale,1,0);
    if(x > *up) {
        return (*up - x) + x;
    }else{
        return x;
    }
}

double F77_SUB(ltgamrnd)(double *shape, double *scale, double *low)
{
    double glow, p=0.0, x=0.0;
    glow=pgamma(*low,*shape,*scale,1,0);
    if(glow>=0.9999) {
        x=(*low)-log(1.0-unif_rand())*(*scale);
    }else{
        p=unif_rand()*(1.0-glow)+glow;
        x=qgamma(p,*shape,*scale,1,0);
    }
    if(x < *low) {
        return (*low - x) + x;
    }else{
        return x;
    }
}

double F77_SUB(invgaussrnd)(double *mu, double *lambda)
{
    double y,r1,r2,m=*mu,la=*lambda;
    y=rchisq(1.0);
    r2=m/(2.0*la)*(2.0*la+m*y+sqrt(4.0*la*m*y+R_pow(m,2)*R_pow(y,2)));
    r1=R_pow(m,2)/r2;
    if(unif_rand() < m/(m+r1)){
        return r1;
    }else{
        return r2;
    }
}

double F77_SUB(betarnd)(double *a, double *b)
{
    return rbeta(*a, *b);
}

double F77_SUB(ltlogisrnd)(double *location, double *scale, double *low)
{
    double u,pcut,x;
    if((*scale)==0.0){
        return ((*location)>(*low)) ? (*location) : (*low);
    }

    pcut=plogis(*low,*location,*scale,1,0);
    if(pcut>0.9999)
        return (*low)+0.0001*(*scale);
    u=unif_rand();
    u=pcut+(1.0-pcut)*u;
    x=qlogis(u,*location,*scale,1,0);

    return x;
}

double F77_SUB(rtlogisrnd)(double *location, double *scale, double *up)
{
    double u,pcut,x;
    if ((*scale)==0.0){
        return ((*location)<(*up)) ? (*location) : (*up);
    }

    pcut=plogis(*up,*location,*scale,1,0);
    if (pcut<0.0001)
        return (*up)-0.0001*(*scale);
    u=unif_rand();
    u=u*pcut;
    x=qlogis(u,*location,*scale,1,0);

    return x;
}

/******************************************************************************************************
 Cumulative density function (CDF)
 ******************************************************************************************************/

 double F77_SUB(cdfnorm)(double *x, double *mu, double *sigma, int *lower_tail, int *give_log)
 {
   return pnorm(*x, *mu, *sigma, *lower_tail, *give_log);
 }

 /******************************************************************************************************
  Probability density function (PDF)
  ******************************************************************************************************/

 double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
 {
 	return dnorm(*x, *mu, *sigma, *give_log);
 }

 double F77_SUB(dexpo)(double *x, double *scale, int *give_log)
 {
     return dexp(*x, *scale, *give_log);
 }

 double F77_SUB(dgamm)(double *x, double *shape, double *scale, int *give_log)
 {
 	return dgamma(*x, *shape, *scale, *give_log);
 }

 double F77_SUB(dald)(double *x, double *location, double *scale, double *p, int *give_log)
 {
     double y=*x,m=*location,s=*scale,tau=*p;
     double consts,exponent=0.0,dens;
     consts=(tau*(1.0-tau)/s);
     if(y<m){
         exponent=exp((1.0-tau)*(y-m)/s);
     }else{
         exponent=exp(-(tau)*(y-m)/s);
     }
     dens=consts*exponent;
     if(*give_log){
         dens=log(dens);
     }
     return dens;
 }


/******************************************************************************************************
 Math functions
 ******************************************************************************************************/

double F77_SUB(gammaln)(double *x)
{
	return lgammafn(*x);
}

double F77_SUB(powerxy)(double *x, double *y)
{
  return R_pow(*x, *y);
}

/******************************************************************************************************
 Print on screen
 ******************************************************************************************************/

void F77_SUB(haprint)(int *func, double *acceptr)
{
  Rprintf("function[%d]: pmet = %.4f > 0.6. Increase metm and redo MCMC loop\n", *func, *acceptr);
}


void F77_SUB(laprint)(int *func, double *acceptr)
{
  Rprintf("function[%d]: pmet = %.4f < 0.3. Reduce metm and redo MCMC loop\n", *func, *acceptr);
}


void F77_SUB(aprint)(int *func, double *acceptr)
{
  Rprintf("function[%d]: pmet = %.4f\n", *func, *acceptr);
}


void F77_SUB(sprint)(int *iter, int *tot_iter, double *sec)
{
	Rprintf("MCMC draws %i of %i (CPU time: %.3f s)\n", *iter, *tot_iter, *sec);
}


/******************************************************************************************************
 Missing value
 ******************************************************************************************************/

int F77_SUB(ismiss)(double*x)
{
  if(R_IsNaN(*x) || R_IsNA(*x)) return 1;
  else return 0;
}
