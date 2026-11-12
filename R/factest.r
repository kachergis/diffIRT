factest=function(object,start=NULL,se=F, control=list()){
 N=object$N
 nit=object$nit
 rt=object$rt
 x=object$score
 con=object$control
 model=object$model
 A=object$A
 W=object$W
 nq=length(A)
 ai=object$par.log[1:nit]
 vi=object$par.log[(nit+1):(2*nit)]
 ter=object$par.log[(2*nit+1):(3*nit)]
 sd2A=object$par.log[3*nit+1]
 sd2V=object$par.log[3*nit+2]
 if(!is.null(start) & length(start)!=(2*N))
   stop("If you provide starting values, this vector should be of size 2 x [number of subjects].\n")
 con=list(eps=if(is.null(object$control$eps)) .01 else object$control$eps,tol=1e-8,halfwidth=NULL,adaptive=TRUE)
 con[names(control)]=control
 if(!is.logical(se)) stop("the 'se' argument should be of type logical")
 
  ou=matrix(,N,2)
  for(i in 1:N){
    xx=calcEZ(mean(x[i,],na.rm=T), var(rt[i,],na.rm=T), mean(rt[i,],na.rm=T), N=ncol(x))
    ou[i,1]=xx[[1]]
    ou[i,2]=xx[[2]]
  }
  ou[is.nan(ou)]=NA
  ou[ou<0]=.1
  ou[,2]=log(ou[,2])
  if(model==2) ou[,1]=log(ou[,1])
   apstart=scale(ou[,2])*sqrt(exp(sd2A))
   vpstart=scale(ou[,1])*sqrt(exp(sd2V))

 rt[is.na(rt)]=-999
 rt[is.na(x)]=-999
 x[is.na(x)]=-999

   pars=c(apstart,vpstart)
  # Adaptive grids from the fit (diffIRT(..., control=list(adaptive=TRUE)), the default):
  # integrate over the other latent variable on each person's grid, and start each
  # search at the person's posterior mean instead of the EZ approximation.
  aq=if(isTRUE(con$adaptive)) object$aq else NULL
  if(!is.null(aq)) pars=c(aq$muA,aq$muV)
  if(!is.null(start)) pars[!is.na(start)]=start[!is.na(start)]

  # Each person's score is a separate 1-D problem (the old joint BFGS over all N
  # parameters scaled as N^2). FACscore minimizes each subject's term of the
  # LLfacA / LLfacV objective with Brent's method, in parallel (see src/diffIRT.c).
  fac = function(which, st, hw) {
    st = as.double(st); st[!is.finite(st)] = 0
    r = .C("FACscore", as.double(rt), as.double(x), as.integer(N), as.integer(nit), as.double(ai),
           as.double(vi), as.double(sd2A), as.double(sd2V), as.double(ter), as.integer(which),
           as.integer(model), as.double(A), as.double(W), as.integer(nq), as.double(con$eps),
           st, as.double(hw), as.double(con$tol), as.integer(!is.null(aq)),
           as.double(if(is.null(aq)) double(N) else if(which==1) aq$muV else aq$muA),
           as.double(if(is.null(aq)) double(N) else if(which==1) aq$sV else aq$sA),
           est=double(N), hess=double(N), conv=integer(N))
    if (any(r$conv != 0))
      warning(sum(r$conv != 0), " subject(s) have a ", c("gamma", "theta")[which],
              " estimate at the edge of the search interval; consider a larger control$halfwidth.")
    r
  }
  hwA = if (is.null(con$halfwidth)) max(3, 10*sqrt(exp(sd2A))) else con$halfwidth
  hwV = if (is.null(con$halfwidth)) max(3, 10*sqrt(exp(sd2V))) else con$halfwidth
  resA = fac(1, pars[1:N], hwA)
  resV = fac(2, pars[(N+1):(2*N)], hwV)
  parA = exp(resA$est)
  parV = if (model == 2) exp(resV$est) else resV$est

  ret=matrix(c(parA,parV),N,2)
  dimnames(ret)[[2]]=c("gamma[p]","theta[p]")
  if(se) {
    # objective is -2 log posterior, so var = 2 / second derivative; delta method for exp()
    varA = 2/resA$hess
    varV = 2/resV$hess
    seA = parA*sqrt(varA)
    seV = if (model == 2) parV*sqrt(varV) else sqrt(varV)
    ret = cbind(parA, seA, parV, seV)
    dimnames(ret)[[2]] = c("gamma[p]","gamma[p].se","theta[p]","theta[p].se")
  }
  return(ret)
  }
 
  