factest=function(object,start=NULL,se=F, control=list(), type=c("MAP","EAP","PV"), npv=5){
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
 con=list(eps=if(is.null(object$control$eps)) .01 else object$control$eps,tol=1e-8,halfwidth=NULL,adaptive=TRUE,threads=object$control$threads,
          nq.eap=21,pv.draws=200,pv.items=NULL)
 con[names(control)]=control
 if(!is.null(con$threads)){ oldThreads=setThreads(con$threads); on.exit(setThreads(oldThreads),add=TRUE) }
 if(!is.logical(se)) stop("the 'se' argument should be of type logical")
 type=match.arg(type)
 if(type!="MAP") return(factestPost(object,type,se,npv,con))
 
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
 
  

# EAP scores and plausible values (see R/posterior.r)
factestPost=function(object,type,se,npv,con){
  N=object$N; nit=object$nit
  rt=object$rt; x=object$score
  rt[is.na(rt)]=-999; rt[is.na(x)]=-999; x[is.na(x)]=-999
  par=object$par.log
  aq=postCentres(object,par,rt,x)
  nms=c("gamma[p]","theta[p]")
  if(type=="EAP"){
    e=eapScores(object,par,rt,x,aq,con$nq.eap)
    if(any(is.na(e[,1]))) warning(sum(is.na(e[,1]))," subject(s) have a zero likelihood everywhere; their scores are NA.")
    if(se){ colnames(e)=c(nms[1],paste0(nms[1],".se"),nms[2],paste0(nms[2],".se")); return(e) }
    e=e[,c(1,3),drop=FALSE]; colnames(e)=nms
    return(e)
  }
  # plausible values: one posterior draw per person per set; with pv.items, each set first
  # draws the item and population parameters from their asymptotic normal distribution
  if(!is.numeric(npv) || length(npv)!=1 || npv<1) stop("'npv' should be a positive integer")
  items=con$pv.items
  if(is.null(items)) items=!is.null(object$hessian)
  if(items && is.null(object$hessian)) stop("pv.items = TRUE needs the covariance matrix of the estimates: re-fit with se = TRUE.")
  if(items){
    V=vcov(object); V=(V+t(V))/2; free=par[match(1:object$npars,object$constrain)]
    ev=eigen(V,symmetric=TRUE)
    if(any(ev$values<=0)) warning("The covariance matrix of the estimates is not positive definite; ",
      "plausible values draw item parameters only along its positive directions.")
    L=ev$vectors%*%diag(sqrt(pmax(ev$values,0)),length(ev$values))
    minRT=log(apply(object$rt,2,min,na.rm=TRUE))   # a draw with Ter above an item's fastest RT
  }                                                # has zero likelihood: draw again
  lowESS=0
  out=lapply(1:npv,function(m){
    p=par; a=aq
    if(items){
      for(tries in 1:100){
        f=free+as.vector(L%*%stats::rnorm(length(free)))
        p[object$constrain!=0]=f[object$constrain[object$constrain!=0]]
        if(all(p[(2*nit+1):(3*nit)]<minRT)) break
      }
      if(tries==100) stop("Could not draw item parameters with every Ter below the item's fastest response time.")
      a=postCentres(object,p,rt,x,fitted=FALSE)
    }
    d=pvDraw(object,p,rt,x,a,con$pv.draws)
    lowESS<<-lowESS+sum(d$ess<10,na.rm=TRUE)
    v=cbind(exp(d$zA),if(object$model==2) exp(d$zV) else d$zV)
    colnames(v)=nms
    v
  })
  if(lowESS>0) warning("For ",lowESS," person-draws the importance sample was small (effective size < 10); consider a larger control$pv.draws.")
  attr(out,"items")=items
  out
}
