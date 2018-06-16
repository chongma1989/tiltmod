#' Empirical power analysis
#' 
#' This function deals with power analysis by calculating the relevant TypeI error,
#' power, and probability of being significant by given global false discovery rate.
#' 
#' @param q A numeric value or a vector of numerical value, represent global false 
#'          discovery rate used for power analysis.
#' @param x A vector of numeric values, represent p-values or left-tail areas 
#'          of test statistics from a differential gene expression study.  
#' @param w A vector of two numeric values, represent the weights  
#'          of the uniform and Beta distributions. See \code{\link[tiltmod]{UBMM}}.
#' @param a A vector of two initial parameter values for Beta distribution. See \code{\link[tiltmod]{UBMM}}.
#' @param precision The precision for convergence. Default value is 1e-8.
#' @param MaxIter The maximum iteration for the EM algorhthm.
#' @param theta A numerical value, represents the exponential tilting parameter 
#'              for the fitted mixture model from x. Defualt is NULL. 
#' @param alpha A numeric value, used to determine the probably null region in method ``m1" 
#'              (see \code{\link[tiltmod]{tqvalue}}). Default is 0.9.
#' @param type A character value, chosen from ``left tail area'' and ``pvalue''. Default is
#'             ``left tail area''.
#' @param rel.tol the accuracy used in \code{\link[stats]{integrate}}. 
#' @param tol the accuracy used in \code{\link[stats]{uniroot}}. 
#' 
#' @return A dataframe consists of q, TypeI, Power, ProbS, respectively. TypeI, Power,
#'         and ProbS are calculated based on the rejection region R(q) and
#'         the empirical mixture model for x. 
#'         \describe{
#'            \item{q}{The global false discovery rates provided in arguments.}
#'            \item{TypeI}{\eqn{P(R(q)|H_0)}}
#'            \item{power}{\eqn{P(R(q)|H_1)}}
#'            \item{ProbS}{\eqn{P(R(q))}}
#'         }
#'         
#'         If theta is provided, then the results contain two data frames as above, one 
#'         is calculated from the non-exponential tilted mixture model and the other 
#'         from the exponential tilted mixture model, respectively.  
#'         
#' @examples
#' x=c(rbeta(50,0.5,0.5),runif(950))
#' q=seq(0.05,0.95,0.05)
#' \dontrun{
#' power(q,x,alpha=0.9,type="left tail area")
#' power(q,x,theta=2,alpha=0.9,type="left tail area")
#' }
#'
#' @export
power<-function(q,x,w=NULL,a=NULL,
                precision=1e-8,MaxIter=10000L,
                theta=NULL,alpha=0.9,
                type=c("left tail area","pvalue"),
                rel.tol=.Machine$double.eps^0.25,
                tol=.Machine$double.eps^0.5){
  n=length(q)
  TypeI=Power=ProbS=rep(0,n)
  eps=.Machine$double.eps^0.75
  
  # fit the mixture model for x
  if(is.null(w) && !is.null(a)) 
    emnull = UBMM(x,a,precision=precision,MaxIter=MaxIter)
  if(!is.null(w) && is.null(a)) 
    emnull = UBMM(x,w,precision=precision,MaxIter=MaxIter)
  if(is.null(w) && is.null(a)) 
    emnull = UBMM(x,precision=precision,MaxIter=MaxIter)
  
  tau = emnull$Weight
  parm = emnull$BetaPar
  
  if(any(is.na(tau) || is.na(parm))) 
    stop("EM algorithm for the mixture model does not converge! 
         May need reinitialize the weights or Beta parameters!")
  
  ## The cdf for the mixture of uniform and Beta distributions
  f <- function(x) tau[1]+tau[2]*dbeta(x,parm[1],parm[2])
  Fw <- function(x) tau[1]*x+tau[2]*pbeta(x,parm[1],parm[2])
  
  ## adjust the mixture model
  if(type %in% "left tail area"){
    lcut=0.5-alpha/2
    rcut=0.5+alpha/2
    A=max(dunif(lcut),dunif(rcut))
    
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
    
    for(i in 1:n){
      lc=tryCatch({
        uniroot(function(x) new.tau[1]*integrate(new.f0,0,x,rel.tol = rel.tol)$value/Fw(x)-q[i],c(eps,0.5),tol = tol)$root
      },error=function(e) NA)
      
      rc=tryCatch({
        uniroot(function(x) new.tau[1]*integrate(new.f0,x,1,rel.tol = rel.tol)$value/(1-Fw(x))-q[i],c(0.5,1-eps),tol = tol)$root
      },error=function(e) NA)
      
      if(is.na(lc) || is.na(rc)){
        TypeI[i]=Power[i]=ProbS[i]=NA
      }else{
        TypeI[i]=integrate(new.f0,0,lc,rel.tol = rel.tol)$value+integrate(new.f0,rc,1,rel.tol = rel.tol)$value
        Power[i]=integrate(new.f1,0,lc,rel.tol = rel.tol)$value+integrate(new.f1,rc,1,rel.tol = rel.tol)$value
        ProbS[i]=Fw(lc)+1-Fw(rc) 
      }
    }
  }
  
  if(type %in% "pvalue"){
    ## adjusted uniform-beta mixture model
    lcut=1-alpha
    A=dunif(lcut)
    
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
    
    for(i in 1:n){
      lc=tryCatch({
        uniroot(function(x) new.tau[1]*integrate(new.f0,0,x,rel.tol = rel.tol)$value/Fw(x)-q[i],c(eps,1-eps),tol = tol)$root
      },error=function(e) NA)
      
      if(is.na(lc)){
        TypeI[i]=Power[i]=ProbS[i]=NA
      }else{
        TypeI[i]=integrate(new.f0,0,lc,rel.tol = rel.tol)$value
        Power[i]=integrate(new.f1,0,lc,rel.tol = rel.tol)$value
        ProbS[i]=Fw(lc)
      }
    }
  }
  
  if(is.null(theta) || !is.numeric(theta)){
    return(data.frame(cbind(q,TypeI,Power,ProbS)))
    if(!is.numeric(theta)) warning("theta must be a numerical value!")
  } 
  
  if(!is.null(theta) && is.numeric(theta)){
    # initialize exponential tilted power 
    tTypeI=tPower=tProbS=rep(0,n)
    h<-function(x) tau[1]/f(x)
    tilt_tau = vector(mode="numeric",2)
    
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
    
    if(type %in% "left tail area"){
      ## adjusted tilted mixture model
      tilt_lcut=uniroot(function(x) tilt_F0(x)-0.5+alpha/2,c(eps,0.5),tol = tol)$root
      tilt_rcut=uniroot(function(x) tilt_F0(x)-0.5-alpha/2,c(0.5,1-eps),tol = tol)$root
      A_t=max(tilt_f0(tilt_lcut),tilt_f0(tilt_rcut))
      
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
      
      for(i in 1:n){
        tlc=tryCatch({
          uniroot(function(x) new.tilt_tau[1]*integrate(new.tilt_f0,0,x,rel.tol = rel.tol)$value/
                    integrate(new.tilt_f,0,x,rel.tol = rel.tol)$value-q[i],c(eps,0.5),tol = tol)$root
        },error=function(e) NA)
        
        trc=tryCatch({
          uniroot(function(x) new.tilt_tau[1]*integrate(new.tilt_f0,x,1,rel.tol = rel.tol)$value/
                    integrate(new.tilt_f,x,1,rel.tol = rel.tol)$value-q[i],c(0.5,1-eps),tol = tol)$root
        },error=function(e) NA)
        
        if(is.na(tlc) || is.na(trc)){
          tTypeI[i]=tPower[i]=tProbS[i]=NA
        }else{
          tTypeI[i]=integrate(new.tilt_f0,0,tlc,rel.tol = rel.tol)$value+integrate(new.tilt_f0,trc,1,rel.tol = rel.tol)$value
          tPower[i]=integrate(new.tilt_f1,0,tlc,rel.tol = rel.tol)$value+integrate(new.tilt_f1,trc,1,rel.tol = rel.tol)$value
          tProbS[i]=integrate(new.tilt_f,0,tlc,rel.tol = rel.tol)$value+integrate(new.tilt_f,trc,1,rel.tol = rel.tol)$value 
        }
      }
    }
    
    if(type %in% "pvalue"){
      ## adjusted tilted mixture model
      tilt_lcut=uniroot(function(x) tilt_F0(x)-1+alpha,c(eps,1),tol = tol)$root
      A_t=tilt_f0(tilt_lcut)
      
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
      
      for(i in 1:n){
        tlc=tryCatch({
          uniroot(function(x) new.tilt_tau[1]*integrate(new.tilt_f0,0,x,rel.tol = rel.tol)$value/
                    integrate(new.tilt_f,0,x,rel.tol = rel.tol)$value-q[i],c(eps,1-eps),tol = tol)$root
        },error=function(e) NA)
     
        if(is.na(tlc) || is.na(trc)){
          tTypeI[i]=tPower[i]=tProbS[i]=NA
        }else{
          tTypeI[i]=integrate(new.tilt_f0,0,tlc,rel.tol = rel.tol)$value
          tPower[i]=integrate(new.tilt_f1,0,tlc,rel.tol = rel.tol)$value
          tProbS[i]=integrate(new.tilt_f,0,tlc,rel.tol = rel.tol)$value
        }
      }
    }
  }
  return(list(data.frame(cbind(q,TypeI,Power,ProbS)),
              data.frame(cbind(q,tTypeI,tPower,tProbS))))
}
