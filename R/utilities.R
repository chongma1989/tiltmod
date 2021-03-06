# Copyright (c) 2016 - 2019 Chong Ma
#
# This file contains the utility functions for the R package. The R package
# is written for reproducible screening of the large scale testing problem
# by combining the cross-validation subsampling and false discovery rate.


#' The exponential tilting function
#'
#' This function tilts the mixture model fitted from the training tail-areas (or p-values)
#' by conditioning on the average of local fdr's from the testing tail-areas (or p-values)
#' @param xl The training left-tail areas (or p-values)
#' @param xt The testing left-tail areas (or p-values)
#' @param f The objective function is tilted. If either xl or f is NULL, f is fitted by
#'    \code{\link[tiltmod]{UBMM}}. Default is NULL.
#' @param h The conditioning function. By default, h = (1-p)/f(x), where 
#'          \eqn{f(x)=(1-p) \times duniform(x)+p \times dbeta(x,\alpha,\beta)}.
#' @param m The constant is used to find the optimal theta such that E(h(x))=m. If m is
#'    NULL, m = mean(h(xt)). Default is NULL.
#' @param interval The interval is used to search the optimal theta. Default is (-100L,100L).
#' @param rel.tol the accuracy used in \code{\link[stats]{integrate}}. 
#' @param ... Arguments to be passed to \code{\link[stats]{uniroot}}.
#' 
#' @return A list includes theta, tau, tilt_tau, tilt_f, tilt_f0, tilt_f1, respectively.
#' @examples
#' xl=c(rbeta(100,0.5,0.5),runif(900))
#' xt=c(rbeta(300,2,3),runif(700))
#'
#' \dontrun{
#' etilt(xl,xt)
#' }
#'
#' @export
etilt <- function(xl,xt,f=NULL,h=NULL,m=NULL,
                  interval=NULL,rel.tol=.Machine$double.eps^0.25,...){
  # fit the mixture model for xl
  if(is.null(f)){
    emnull = UBMM(xl)
    tau = emnull$Weight
    parm = emnull$BetaPar
    f <- function(x) tau[1]+tau[2]*dbeta(x,parm[1],parm[2])
    tilt_tau = vector(mode="numeric",2)
  }

  # entropy function
  if(is.null(h)) h <- function(x) tau[1]/f(x)

  # constant
  if(is.null(m)) m = mean(h(xt))

  # calculate the optimal theta such that E[h(x)]=m
  TiltGmean<-function(theta){
    numerator <- function(x,the) h(x)*exp(the*h(x))*f(x)
    denomenator <- function(x,the) exp(the*h(x))*f(x)

    return(sapply(theta,function(t)
      log(integrate(numerator,0,1,t,rel.tol = rel.tol,stop.on.error = FALSE)$value)-
        log(integrate(denomenator,0,1,t,rel.tol = rel.tol,stop.on.error=FALSE)$value)-log(m)))
  }

  # calculate the tilted distribution
  if(is.null(interval)) interval = c(-100L,100L)
  theta=tryCatch(uniroot(TiltGmean,interval=interval,...)$root,
                 warning=function(w) {
                   warning("The optimal theta is not found! Theta is set 0 without tilting!");
                   return(0)
                   })

  tilt_kf <- function(x) exp(theta*h(x))*f(x)
  tilt_kf0 <- function(x) exp(theta*h(x))*dunif(x)
  tilt_kf1 <- function(x) exp(theta*h(x))*dbeta(x,parm[1],parm[2])

  cnst = integrate(tilt_kf,0,1,stop.on.error = FALSE)$value
  cnst0 = integrate(tilt_kf0,0,1,stop.on.error = FALSE)$value
  cnst1 = integrate(tilt_kf1,0,1,stop.on.error = FALSE)$value

  tilt_f <- function(x) tilt_kf(x)/cnst
  tilt_f0 <- function(x) tilt_kf0(x)/cnst0
  tilt_f1 <- function(x) tilt_kf1(x)/cnst1

  tilt_tau[1] = tau[1]*cnst0/cnst
  tilt_tau[2] = 1 - tilt_tau[1]

  return(list(theta = theta,
              tau = tau,
              tilt_tau = tilt_tau,
              tilt_f = tilt_f,
              tilt_f0 = tilt_f0,
              tilt_f1 = tilt_f1))
}

#' The exponential tilting mixture model
#'
#' This function tilts the mixture model fitted from the training tail-areas (or p-values)
#' by conditioning on the average of local fdr's from the testing tail-areas (or p-values)
#' @param xl The training left-tail areas (or p-values)
#' @param xt The testing left-tail areas (or p-values)
#' @param w A vector of two numeric values, representing the weights  
#'          of the uniform and Beta distributions. See \code{\link[tiltmod]{UBMM}}.
#' @param a A vector of two initial parameter values for Beta distribution. See \code{\link[tiltmod]{UBMM}}.
#' @param precision The precision for convergence. Default value is 1e-8.
#' @param MaxIter The maximum iteration for the EM algorhthm.
#' @param interval A vector of two numeric values, which determines the range to
#'                 search the optimal theta. Default is c(-1000L,1000L).
#' @param adjust Whether or not to do the model adjustment. Default is TRUE.
#' @param method A character chosen from m1, m2. Default is m1.
#' @param type A character value, chosen from ``left tail area'' and ``pvalue''. Default is
#'             ``left tail area''.
#' @param alpha A numeric value. Used in method ``m1'' to determine the 
#'              probably null region. Default is 0.9.
#' @param q A numeric value. The global false discovery rate used in method ``m2'', 
#'          to determine the probable null region. Default is 0.1.
#' @param ncores The number of cpus used for implementing this function.
#' @param rel.tol the accuracy used in \code{\link[stats]{integrate}}. 
#' @param tol the accuracy used in \code{\link[stats]{uniroot}}. 
#' @param eps the smallest positive precision. If x < eps, x = eps; if x > 1-eps, x = 1-eps.
#' 
#' @return A dataframe includes xl, xt, fdr, FDR, tfdr, and tFDR, respectively. fdr and
#'         FDR are the local and global false discovery rate for each value of xt.
#'         tfdr and tFDR are the corresponding tilted local and global false discovery
#'         rate, respectively.
#'
#'         The optimal theta calculated by solving log(E(exp(thetah(x))))-ctheta,
#'         where c=mean(h(xt)).
#' @examples
#' xl=c(rbeta(50,0.2,0.2),runif(950))
#' xt=c(rbeta(50,0.1,0.1),runif(950))
#'
#' \dontrun{
#' tqvalue(xl,xt,ncores=4,adjust=FALSE,type="left tail area")
#' }
#'
#' @export
tqvalue <- function(xl,xt,
                    w=NULL,a=NULL,
                    precision=1e-8,MaxIter=10000L,
                    interval=NULL,adjust=TRUE,
                    method=c("m1","m2"),
                    type=c("left tail area","pvalue"),
                    alpha=0.9,q=0.1,ncores=1,
                    rel.tol=.Machine$double.eps^0.25,
                    tol=.Machine$double.eps^0.5,
                    eps=.Machine$double.eps^0.6){
  # fit the mixture model for xl
  if(is.null(w) && !is.null(a)) 
    emnull = UBMM(xl,a,precision=precision,MaxIter=MaxIter)
  if(!is.null(w) && is.null(a)) 
    emnull = UBMM(xl,w,precision=precision,MaxIter=MaxIter)
  if(!is.null(w) && !is.null(a)) 
    emnull = UBMM(xl,w,a,precision=precision,MaxIter=MaxIter)
  if(is.null(w) && is.null(a)) 
    emnull = UBMM(xl,precision=precision,MaxIter=MaxIter)
  
  tau = emnull$Weight
  parm = emnull$BetaPar
  tilt_tau = vector(mode="numeric",2)
  
  if(any(is.na(tau) || is.na(parm))) 
    stop("EM algorithm for the mixture model does not converge! 
       May need reinitialize the weights or Beta parameters!")
  
  ## The pdf and cdf for the mixture of uniform and Beta distributions
  f <- function(x) tau[1]+tau[2]*dbeta(x,parm[1],parm[2])
  Fw <- function(x) tau[1]*x+tau[2]*pbeta(x,parm[1],parm[2])
  
  ## Tilting function
  h<-function(x) tau[1]/f(x)
  m=mean(h(xt))
  
  # calculate the optimal theta such that E[h(x)]=m
  TiltGmean<-function(theta){
    numerator <- function(x,the) h(x)*exp(the*h(x))*f(x)
    denomenator <- function(x,the) exp(the*h(x))*f(x)
    
    return(sapply(theta,function(t)
      log(integrate(numerator,0,1,t,rel.tol = rel.tol,stop.on.error = FALSE)$value)-
        log(integrate(denomenator,0,1,t,rel.tol = rel.tol,stop.on.error=FALSE)$value)-log(m)))
  }
  
  # calculate the tilted distribution
  # check interval is appropriate
  if(is.null(interval)) interval = c(-1000L,1000L)
  while(interval[1] < -100L || interval[2] > 100L){
    interval=tryCatch({
      if(prod(TiltGmean(interval))<0){
        break
      }else{
        c(interval[1]+100L,interval[2]-100L)
      }}, error=function(e) c(interval[1]+100L,interval[2]-100L))
  }
  
  theta=tryCatch(uniroot(TiltGmean,interval=interval)$root,
                 warning=function(w) {
                   warning("The optimal theta is not found! Theta is set 0 without tilting!");
                   return(0)
                 })
  
  tilt_kf <- function(x) exp(theta*h(x))*f(x)
  tilt_kf0 <- function(x) exp(theta*h(x))*dunif(x)
  tilt_kf1 <- function(x) exp(theta*h(x))*dbeta(x,parm[1],parm[2])
  
  cnst = integrate(tilt_kf,0,1,stop.on.error = FALSE)$value
  cnst0 = integrate(tilt_kf0,0,1,stop.on.error = FALSE)$value
  cnst1 = integrate(tilt_kf1,0,1,stop.on.error = FALSE)$value
  
  tilt_f <- function(x) tilt_kf(x)/cnst
  tilt_f0 <- function(x) tilt_kf0(x)/cnst0
  tilt_f1 <- function(x) tilt_kf1(x)/cnst1
  
  tilt_tau[1] = tau[1]*cnst0/cnst
  tilt_tau[2] = 1 - tilt_tau[1]
  
  tilt_F0 <- function(x) integrate(tilt_f0,0,x,rel.tol = rel.tol,stop.on.error = FALSE)$value
  tilt_F <- function(x) integrate(tilt_f,0,x,rel.tol = rel.tol,stop.on.error = FALSE)$value
  
  ## use untilted distribution to find the symetric piont t to x
  froot1 <- function(t,x) Fw(t)+Fw(x)-1
  ## use tilted distribution to find the symetric piont t to x
  froot2 <- function(t,x) {
    return(ifelse(t!=0,integrate(tilt_kf,0,t,rel.tol = rel.tol,stop.on.error = FALSE)$value,0)-
             ifelse(x!=1,integrate(tilt_kf,x,1,rel.tol = rel.tol,stop.on.error = FALSE)$value,0))
  }
  
  ## error checking
  if(all(method %in% c("m1","m2")) && length(method)==2) method="m1"
  if(all(type %in% c("left tail area","pvalue")) && length(type)==2) type="left tail area"
  if(length(method) > 1 || !(method %in% c("m1","m2"))) stop("Must choose one method from m1 or m2!")
  if(length(type) > 1 || !(type %in% c("left tail area","pvalue"))) stop("Must choose one type from left tail area and pvalue!")
  
  if(type %in% "left tail area"){
    if(method %in% "m1"){
      ## adjusted uniform-beta mixture model
      lcut=0.5-alpha/2
      rcut=0.5+alpha/2
      A=max(dunif(lcut),dunif(rcut))
      
      ## adjusted tilted mixture model
      tilt_lcut=uniroot(function(x) tilt_F0(x)-0.5+alpha/2,c(eps,0.5),tol = tol)$root
      tilt_rcut=uniroot(function(x) tilt_F0(x)-0.5-alpha/2,c(0.5,1-eps),tol = tol)$root
      A_t=max(tilt_f0(tilt_lcut),tilt_f0(tilt_rcut))
    }
    
    if(method %in% "m2"){
      ## adjusted uniform-beta mixture model
      tryCatch({
        lcut=uniroot(function(x) tau[1]*x/Fw(x)-q,c(eps,0.5),tol = tol)$root
        rcut=uniroot(function(x) tau[1]*(1-x)/(1-Fw(x))-q,c(0.5,1-eps),tol = tol)$root
      },error=function(e) stop("The method 2 is not applicable! Try method 1!"))
      A=max(dunif(lcut),dunif(rcut))
      
      #adjusted tilted mixture model
      tryCatch({
        tilt_lcut=uniroot(function(x) tilt_tau[1]*tilt_F0(x)/tilt_F(x)-q,c(eps,0.5),tol = tol)$root
        tilt_rcut=uniroot(function(x) tilt_tau[1]*integrate(tilt_f0,x,1)$value/integrate(tilt_f,x,1)$value-q,c(0.5,1-eps),tol = tol)$root
      },error=function(e) stop("The method 2 is not applicable! Try method 1!"))
      A_t=max(tilt_f0(tilt_lcut),tilt_f0(tilt_rcut))
    }
    
    ## adjusted uniform-beta mixture model
    f11<-function(x) ifelse(x>=lcut,ifelse(x<=rcut,0,dbeta(x,parm[1],parm[2])-A),
                           dbeta(x,parm[1],parm[2])-A)
    f10<-function(x) ifelse(x>=lcut,ifelse(x<=rcut,dbeta(x,parm[1],parm[2]),A),A)
    
    f1<-function(x) ifelse(f11(x)<0,0,f11(x))
    f0<-function(x) tau[1]*dunif(x)+tau[2]*f10(x)
    
    A1=integrate(f1,0,1,rel.tol = rel.tol)$value
    A0=integrate(f0,0,1,rel.tol = rel.tol)$value
    new.tau=c(1-tau[2]*A1,tau[2]*A1)
    
    new.f1<-function(x) f1(x)/A1
    new.f0<-function(x) f0(x)/A0
    new.f=f
    
    #adjusted tilted mixture model
    tilt_f111<-function(x) ifelse(x>=tilt_lcut,ifelse(x<=tilt_rcut,0,tilt_f1(x)-A_t),
                                 tilt_f1(x)-A_t)
    tilt_f110<-function(x) ifelse(x>=tilt_lcut,ifelse(x<=tilt_rcut,tilt_f1(x),A_t),A_t)
    
    tilt_f11<-function(x) ifelse(tilt_f111(x)<0,0,tilt_f111(x))
    tilt_f10<-function(x) tilt_tau[1]*tilt_f0(x)+tilt_tau[2]*tilt_f110(x)
    
    tilt_A1=integrate(tilt_f11,0,1,rel.tol = rel.tol)$value
    tilt_A0=integrate(tilt_f10,0,1,rel.tol = rel.tol)$value
    new.tilt_tau=c(1-tilt_tau[2]*tilt_A1,tilt_tau[2]*tilt_A1)
    
    new.tilt_f1<-function(x) tilt_f11(x)/tilt_A1
    new.tilt_f0<-function(x) tilt_f10(x)/tilt_A0
    new.tilt_f=tilt_f
  }
  
  if(type %in% "pvalue"){
    if(method %in% "m1"){
      ## adjusted uniform-beta mixture model
      lcut=1-alpha
      A=dunif(lcut)
      
      ## adjusted tilted mixture model
      tilt_lcut=uniroot(function(x) tilt_F0(x)-1+alpha,c(eps,1),tol = tol)$root
      A_t=tilt_f0(tilt_lcut)
    }
    
    if(method %in% "m2"){
      ## adjusted uniform-beta mixture model
      tryCatch({
        lcut=uniroot(function(x) tau[1]*x/Fw(x)-q,c(eps,0.5),tol = tol)$root
      },error=function(e) stop("The method 2 is not applicable! Try method 1!"))
      A=dunif(lcut)
      
      #adjusted tilted mixture model
      tryCatch({
        tilt_lcut=uniroot(function(x) tilt_tau[1]*tilt_F0(x)/tilt_F(x)-q,c(eps,0.5),tol = tol)$root
      },error=function(e) stop("The method 2 is not applicable! Try method 1!"))
      A_t=tilt_f0(tilt_lcut)
    }
    
    ## adjusted uniform-beta mixture model
    f11<-function(x) ifelse(x>lcut,0,dbeta(x,parm[1],parm[2])-A)
    f10<-function(x) ifelse(x>lcut,dbeta(x,parm[1],parm[2]),A)
    
    f1<-function(x) ifelse(f11(x)<0,0,f11(x))
    f0<-function(x) tau[1]*dunif(x)+tau[2]*f10(x)
    
    A1=integrate(f1,0,1,rel.tol = rel.tol)$value
    A0=integrate(f0,0,1,rel.tol = rel.tol)$value
    new.tau=c(1-tau[2]*A1,tau[2]*A1)
    
    new.f1<-function(x) f1(x)/A1
    new.f0<-function(x) f0(x)/A0
    new.f=f
    
    #adjusted tilted mixture model
    tilt_f111<-function(x) ifelse(x>tilt_lcut,0,tilt_f1(x)-A_t)
    tilt_f110<-function(x) ifelse(x>tilt_lcut,tilt_f1(x),A_t)
    
    tilt_f11<-function(x) ifelse(tilt_f111(x)<0,0,tilt_f111(x))
    tilt_f10<-function(x) tilt_tau[1]*tilt_f0(x)+tilt_tau[2]*tilt_f110(x)
    
    tilt_A1=integrate(tilt_f11,0,1,rel.tol = rel.tol)$value
    tilt_A0=integrate(tilt_f10,0,1,rel.tol = rel.tol)$value
    new.tilt_tau=c(1-tilt_tau[2]*tilt_A1,tilt_tau[2]*tilt_A1)
    
    new.tilt_f1<-function(x) tilt_f11(x)/tilt_A1
    new.tilt_f0<-function(x) tilt_f10(x)/tilt_A0
    new.tilt_f=tilt_f
  }
  
  i=NULL ## initializer for iterations
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  if(type %in% "left tail area"){
    # using left tail area
    if(adjust){
      ttTab <- foreach::foreach(i=1:length(xt),.combine=rbind,.multicombine=T,
                       .packages = c('foreach','stats')) %dopar%{
                         tryCatch({
                           ## Non-tilting method
                           ## calculate the reflex point to the observed ptvalue
                           sym.point=tryCatch(uniroot(froot1,xt[i],interval=c(0,1),tol=tol)$root,
                                              error=function(e) 1-xt[i])
                           
                           left=min(xt[i],sym.point)  ##denote the left point wrt "pvalue"
                           right=max(xt[i],sym.point) ##denote the right point wrt "pvalue"
                           
                           ##probability of rejection based on untilt.f0 and untilt.fw
                           if(left < eps && right > 1-eps){
                             F0_x=tryCatch({
                               ifelse(is.na(new.f0(xt[i])),new.f0(eps)+new.f0(1-eps),new.f0(xt[i]))
                             },error=function(e) new.f0(eps)+new.f0(1-eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(new.f(xt[i])),new.f(eps)+new.f(1-eps),new.f(xt[i]))
                             },error=function(e) new.f(eps)+new.f(1-eps))
                           }else if(left < eps && right < 1-eps){
                             F0_x=2*integrate(new.f0,right,1,stop.on.error = FALSE)$value
                             Fw_x=2*integrate(new.f,right,1,stop.on.error = FALSE)$value
                           }else if(left > eps && right > 1-eps){
                             F0_x=2*integrate(new.f0,0,left,stop.on.error = FALSE)$value
                             Fw_x=2*integrate(new.f,0,left,stop.on.error = FALSE)$value
                           }else {
                             F0_x=integrate(new.f0,0,left,stop.on.error = FALSE)$value+
                               integrate(new.f0,right,1,stop.on.error = FALSE)$value
                             Fw_x=integrate(new.f,0,left,stop.on.error = FALSE)$value+
                               integrate(new.f,right,1,stop.on.error = FALSE)$value
                           }
                           
                           fdr=new.tau[1]*new.f0(xt[i])/new.f(xt[i])
                           FDR=new.tau[1]*F0_x/Fw_x
                           
                           ## Tilting method
                           sym.pointt=tryCatch(uniroot(froot2,xt[i],interval=c(0,1),tol=tol)$root,
                                               error=function(e) 1-xt[i])
                           tleft=min(xt[i],sym.pointt)  ##denote the left point wrt "pvalue"
                           tright=max(xt[i],sym.pointt) ##denote the right point wrt "pvalue"
                           
                           ##probability of rejection based on tilt.f0 and tilt.fw
                           if(tleft < eps && tright > 1-eps){
                             tF0_x=tryCatch({
                               ifelse(is.na(new.tilt_f0(xt[i])),new.tilt_f0(eps)+new.tilt_f0(1-eps),new.tilt_f0(xt[i]))
                             },error=function(e) new.tilt_f0(eps)+new.tilt_f0(1-eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(new.tilt_f(xt[i])),new.tilt_f(eps)+new.tilt_f(1-eps),new.tilt_f(xt[i]))
                             },error=function(e) new.tilt_f(eps)+new.tilt_f(1-eps))
                           }else if(tleft < eps && tright < 1-eps){
                             tF0_x=2*integrate(new.tilt_f0,tright,1,stop.on.error = FALSE)$value
                             tFw_x=2*integrate(new.tilt_f,tright,1,stop.on.error = FALSE)$value
                           }else if(tleft > eps && tright > 1-eps){
                             tF0_x=2*integrate(new.tilt_f0,0,tleft,stop.on.error = FALSE)$value
                             tFw_x=2*integrate(new.tilt_f,0,tleft,stop.on.error = FALSE)$value
                           }else {
                             tF0_x=integrate(new.tilt_f0,0,tleft,stop.on.error = FALSE)$value+
                               integrate(new.tilt_f0,tright,1,stop.on.error = FALSE)$value
                             tFw_x=integrate(new.tilt_f,0,tleft,stop.on.error = FALSE)$value+
                               integrate(new.tilt_f,tright,1,stop.on.error = FALSE)$value
                           }
                           
                           tfdr=new.tilt_tau[1]*new.tilt_f0(xt[i])/new.tilt_f(xt[i])
                           tFDR=new.tilt_tau[1]*tF0_x/tFw_x
                           
                           c(fdr=min(abs(fdr),1),FDR=min(abs(FDR),1),
                             tfdr=min(abs(tfdr),1),tFDR=min((tFDR),1))
                         },
                         error=function(e) NA
                         )
                       }
    }else{
      ttTab <- foreach::foreach(i=1:length(xt),.combine=rbind,.multicombine=T,
                       .packages = c('foreach','stats')) %dopar%{
                         tryCatch({
                           ## Non-tilting method
                           ## calculate the reflex point to the observed ptvalue
                           sym.point=tryCatch(uniroot(froot1,xt[i],interval=c(0,1),tol=tol)$root,
                                              error=function(e) 1-xt[i])
                           left=min(xt[i],sym.point)  ##denote the left point wrt "pvalue"
                           right=max(xt[i],sym.point) ##denote the right point wrt "pvalue"
                           
                           ##probability of rejection based on untilt.f0 and untilt.fw
                           if(left < eps && right > 1-eps){
                             F0_x=tryCatch({
                               ifelse(is.na(dunif(xt[i])),dunif(eps)+dunif(1-eps),dunif(xt[i]))
                             },error=function(e) dunif(eps)+dunif(1-eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(f(xt[i])),f(eps)+f(1-eps),f(xt[i]))
                             },error=function(e) f(eps)+f(1-eps))
                           }else if(left < eps && right < 1-eps){
                             F0_x=2*(1-right)
                             Fw_x=2*(1-Fw(right))
                           }else if(left > eps && right > 1-eps){
                             F0_x=2*left
                             Fw_x=2*Fw(left)
                           }else {
                             F0_x=left+1-right
                             Fw_x=Fw(left)+1-Fw(right)
                           }
                           
                           fdr=tau[1]/f(xt[i])
                           FDR=tau[1]*F0_x/Fw_x
                           
                           ## Tilting method
                           sym.pointt=tryCatch(uniroot(froot2,xt[i],interval=c(0,1),tol=tol)$root,
                                               error=function(e) 1-xt[i])
                           tleft=min(xt[i],sym.pointt)  ##denote the left point wrt "pvalue"
                           tright=max(xt[i],sym.pointt) ##denote the right point wrt "pvalue"
                           
                           ##probability of rejection based on tilt.f0 and tilt.f
                           if(tleft < eps && tright > 1-eps){
                             tF0_x=tryCatch({
                               ifelse(is.na(tilt_f0(xt[i])),tilt_f0(eps)+tilt_f0(1-eps),tilt_f0(xt[i]))
                             },error=function(e) tilt_f0(eps)+tilt_f0(1-eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(tilt_f(xt[i])),tilt_f(eps)+tilt_f(1-eps),tilt_f(xt[i]))
                             },error=function(e) tilt_f(eps)+tilt_f(1-eps))
                           }else if(tleft < eps && tright < 1-eps){
                             tF0_x=2*integrate(tilt_f0,tright,1,stop.on.error = FALSE)$value
                             tFw_x=2*integrate(tilt_f,tright,1,stop.on.error = FALSE)$value
                           }else if(tleft > eps && tright > 1-eps){
                             tF0_x=2*integrate(tilt_f0,0,tleft,stop.on.error = FALSE)$value
                             tFw_x=2*integrate(tilt_f,0,tleft,stop.on.error = FALSE)$value
                           }else {
                             tF0_x=integrate(tilt_f0,0,tleft,stop.on.error = FALSE)$value+
                               integrate(tilt_f0,tright,1,stop.on.error = FALSE)$value
                             tFw_x=integrate(tilt_f,0,tleft,stop.on.error = FALSE)$value+
                               integrate(tilt_f,tright,1,stop.on.error = FALSE)$value
                           }
                           
                           tfdr=tilt_tau[1]*tilt_f0(xt[i])/tilt_f(xt[i])
                           tFDR=tilt_tau[1]*tF0_x/tFw_x
                           
                           c(fdr=min(abs(fdr),1),FDR=min(abs(FDR),1),
                             tfdr=min(abs(tfdr),1),tFDR=min((tFDR),1))
                         },
                         error=function(e) NA
                         )
                       }
    }
  }
  
  if(type %in% "pvalue"){
    # using pvalues
    if(adjust){
      ttTab <- foreach::foreach(i=1:length(xt),.combine=rbind,.multicombine=T,
                       .packages = c('foreach','stats')) %dopar%{
                         tryCatch({
                           ## Non-tilting & tilting methods
                           if(xt[i] < eps){
                             F0_x=tryCatch({
                               ifelse(is.na(new.f0(xt[i])),new.f0(eps),new.f0(xt[i]))
                             },error=function(e) new.f0(eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(new.f(xt[i])),new.f(eps),new.f(xt[i]))
                             },error=function(e) new.f(eps))
                             
                             tF0_x=tryCatch({
                               ifelse(is.na(new.tilt_f0(xt[i])),new.tilt_f0(eps),new.tilt_f0(xt[i]))
                             },error=function(e) new.tilt_f0(eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(new.tilt_f(xt[i])),new.tilt_f(eps),new.tilt_f(xt[i]))
                             },error=function(e) new.tilt_f(eps))
                           } else if(xt[i] > 1-eps){
                             F0_x=tryCatch({
                               ifelse(is.na(new.f0(xt[i])),new.f0(1-eps),new.f0(xt[i]))
                             },error=function(e) new.f0(1-eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(new.f(xt[i])),new.f(1-eps),new.f(xt[i]))
                             },error=function(e) new.f(1-eps))
                             
                             tF0_x=tryCatch({
                               ifelse(is.na(new.tilt_f0(xt[i])),new.tilt_f0(1-eps),new.tilt_f0(xt[i]))
                             },error=function(e) new.tilt_f0(1-eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(new.tilt_f(xt[i])),new.tilt_f(1-eps),new.tilt_f(xt[i]))
                             },error=function(e) new.tilt_f(1-eps))
                           } else{
                             F0_x=integrate(new.f0,0,xt[i],stop.on.error = FALSE)$value
                             Fw_x=integrate(new.f,0,xt[i],stop.on.error = FALSE)$value
                             
                             tF0_x=integrate(new.tilt_f0,0,xt[i],stop.on.error = FALSE)$value
                             tFw_x=integrate(new.tilt_f,0,xt[i],stop.on.error = FALSE)$value
                           }
                           
                           fdr=new.tau[1]*new.f0(xt[i])/new.f(xt[i])
                           FDR=new.tau[1]*F0_x/Fw_x
                           
                           tfdr=new.tilt_tau[1]*new.tilt_f0(xt[i])/new.tilt_f(xt[i])
                           tFDR=new.tilt_tau[1]*tF0_x/tFw_x
                           
                           c(fdr=min(abs(fdr),1),FDR=min(abs(FDR),1),
                             tfdr=min(abs(tfdr),1),tFDR=min(abs(tFDR),1))
                         },
                         error=function(e) NA
                         )
                       }
    }else{
      ttTab <- foreach::foreach(i=1:length(xt),.combine=rbind,.multicombine=T,
                       .packages = c('foreach','stats')) %dopar%{
                         tryCatch({
                           ## Non-tilting & tilting methods
                           if(xt[i] < eps){
                             F0_x=tryCatch({
                               ifelse(is.na(dunif(xt[i])),dunif(eps),dunif(xt[i]))
                             },error=function(e) dunif(eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(f(xt[i])),f(eps),f(xt[i]))
                             },error=function(e) f(eps))
                             
                             tF0_x=tryCatch({
                               ifelse(is.na(tilt_f0(xt[i])),tilt_f0(eps),tilt_f0(xt[i]))
                             },error=function(e) tilt_f0(eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(tilt_f(xt[i])),tilt_f(eps),tilt_f(xt[i]))
                             },error=function(e) tilt_f(eps))
                           } else if(xt[i] > 1-eps){
                             F0_x=tryCatch({
                               ifelse(is.na(dunif(xt[i])),dunif(1-eps),dunif(xt[i]))
                             },error=function(e) dunif(1-eps))
                             Fw_x=tryCatch({
                               ifelse(is.na(f(xt[i])),f(1-eps),f(xt[i]))
                             },error=function(e) f(1-eps))
                             
                             tF0_x=tryCatch({
                               ifelse(is.na(tilt_f0(xt[i])),tilt_f0(1-eps),tilt_f0(xt[i]))
                             },error=function(e) tilt_f0(1-eps))
                             tFw_x=tryCatch({
                               ifelse(is.na(tilt_f(xt[i])),tilt_f(1-eps),tilt_f(xt[i]))
                             },error=function(e) tilt_f(1-eps))
                           } else{
                             F0_x=xt[i]
                             Fw_x=Fw(xt[i])
                             
                             tF0_x=integrate(tilt_f0,0,xt[i],stop.on.error = FALSE)$value
                             tFw_x=integrate(tilt_f,0,xt[i],stop.on.error = FALSE)$value
                           }
                          
                           fdr=tau[1]/f(xt[i])
                           FDR=tau[1]*F0_x/Fw_x
                           
                           tfdr=tilt_tau[1]*tilt_f0(xt[i])/tilt_f(xt[i])
                           tFDR=tilt_tau[1]*tF0_x/tFw_x
                           
                           c(fdr=min(abs(fdr),1),FDR=min(abs(FDR),1),
                             tfdr=min(abs(tfdr),1),tFDR=min(abs(tFDR),1))
                         },
                         error=function(e) NA
                         )
                       }
    }
  }
  parallel::stopCluster(cl)
  rownames(ttTab)=NULL
  return(list(ttTable=data.frame(xl=xl,xt=xt,ttTab),theta=theta))
}


#' Frequency network
#'
#' This function displays the frequency network of discovered differentially expressed
#' genes.
#'
#' @param x A list of discovered differentially expressed genes
#' @param Simplify Logical indicating whether to discard the genes with lower relative
#'                 freqency than the threshold. Default is FALSE.
#' @param threshold A numeric value determining the cutoff point, where the genes
#'                  are discarded with lower relative frequency than it. Default is 0.05.
#' @param max.ew A numeric value. The maximum edge width in the network plot.
#' @param max.size A numeric value, displays the maximum node size. 
#' @param node A logical value, illustrates the feature name in a circle. Default is TRUE. 
#' @param directed A logical value, indicates whether the edges are shown in directions. Default 
#'                 is FALSE.
#' @param print.graph A logical value. Default is TRUE. If FALSE, then do not print 
#'                    the frequency network.  
#' @param ... Arguments to be passed to \code{\link[igraph]{plot.igraph}}.
#' 
#' @return A network plot.
#' @import plyr,utils,igraph
#'
#' @examples
#' x=list(c(1,3,4),c(2,4,5),c(3,5,1),c(4,1,2),c(5,2,3))
#' \dontrun{
#' require(igraph)
#' fnet(x)
#' fnet(x,layout=layout.fruchterman.reingold,vertex.color="grey60")
#' }
#'
#' @export
fnet <- function(x,Simplify=FALSE,threshold=0.05,
                 max.ew=3,max.size=24,node=TRUE,
                 directed=FALSE,print.graph=TRUE,...)
  {
  if(!is.list(x)) stop("x should be a list of significant cases (genes)!")
  x=lapply(x, function(t) as.character(t))
  n = length(x)
  xc = plyr::count(unlist(x))
  xc = xc[order(xc$freq,decreasing = TRUE),]
  xc$rfreq = xc$freq/n
  xc$size=max.size-(1-xc$rfreq)^0.4*10
  colnames(xc) = c("name","freq","rfreq","size")
  rownames(xc) = NULL
  
  if(Simplify || threshold !=0){
    unique_name = as.character(with(xc,name[rfreq >= threshold]))
    x = lapply(x, function(t) t[!is.na(match(t,unique_name))])
    xc = xc[match(unique_name,as.character(xc$name)),]
  }
  
  ##create edges
  create_edges<-function(x){
    if(length(x)>1){
      return(c(combn(x,2)))
    }else{
      return(rep(x,2))
    }
  }
  
  x_edges <- unlist(lapply(x, function(t) create_edges(t)))
  x_net <- igraph::graph(x_edges,directed = directed)
  igraph::E(x_net)$weight <- rep(1,length(igraph::E(x_net)))
  
  x_nets <- igraph::simplify(x_net,remove.multiple=TRUE,remove.loops=TRUE,
                     edge.attr.comb=c(weight="sum"))
  igraph::V(x_nets)$size <- xc$size[match(as.character(igraph::V(x_nets)$name),
                                          as.character(xc$name))]
  
  ew = rank(igraph::E(x_nets)$weight)
  ew = 1+ew/length(ew)*max.ew
  if(node){
    vsize=igraph::V(x_nets)$size
  }else{
    vsize=0
  }
  
  fTab=xc[,c("name","freq","rfreq")]
  if(print.graph){
    igraph::plot.igraph(x_nets,
         vertex.size=vsize,
         edge.width=ew,...) 
  }
  
  return(fTab)
}
