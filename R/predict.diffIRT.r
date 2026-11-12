# Expected accuracy and mean response time under a fitted diffIRT model.
#
# For a person with boundary gamma and drift theta on item i: a = gamma/a[i];
# v = theta - v[i] (D model) or theta/v[i] (Q model). With an unbiased start (z = a/2) and
# unit diffusion constant,
#   P(correct) = 1 / (1 + exp(-a v)),
#   E[RT]      = Ter[i] + a/(2v) * tanh(a v / 2)     (a^2/4 + Ter[i] when v = 0).

diffExpect=function(gamma,theta,ai,vi,ter,model){
  a=outer(gamma,ai,"/")
  v=if(model==1) outer(theta,vi,"-") else outer(theta,vi,"/")
  av=a*v
  prob=1/(1+exp(-av))
  mdt=ifelse(abs(v)<1e-8, a^2/4, a/(2*v)*tanh(av/2))
  list(prob=prob, rt=sweep(mdt,2,ter,"+"))
}

predict.diffIRT=function(object,gamma=NULL,theta=NULL,...){
  nit=object$nit
  ai=object$par[1:nit]; vi=object$par[(nit+1):(2*nit)]; ter=object$par[(2*nit+1):(3*nit)]
  oA=object$par[3*nit+1]; oV=object$par[3*nit+2]
  items=paste("item",1:nit)
  if(is.null(gamma) && is.null(theta)){
    # population-level expectations: integrate over the person distribution
    A=object$A; W=object$W
    g=rep(exp(oA*A),length(A)); l=rep(oV*A,each=length(A)); w=as.vector(outer(W,W))
    e=diffExpect(g, if(object$model==1) l else exp(l), ai, vi, ter, object$model)
    return(data.frame(prob=colSums(w*e$prob), rt=colSums(w*e$rt), row.names=items))
  }
  if(is.null(gamma) || is.null(theta) || length(gamma)!=length(theta))
    stop("Provide both 'gamma' and 'theta', of equal length (or neither, for population-level expectations).")
  if(any(gamma<=0)) stop("'gamma' should be positive.")
  if(object$model==2 && any(theta<=0)) stop("In the Q-diffusion model 'theta' should be positive.")
  e=diffExpect(gamma, theta, ai, vi, ter, object$model)
  dimnames(e$prob)=dimnames(e$rt)=list(NULL,items)
  e
}

fitted.diffIRT=function(object,...){
  fs=factest(object,...)
  predict.diffIRT(object, gamma=fs[,"gamma[p]"], theta=fs[,"theta[p]"])
}

simulate.diffIRT=function(object,nsim=1,seed=NULL,keep.missing=TRUE,...){
  if(!is.null(seed)) set.seed(seed)
  nit=object$nit; N=object$N
  ai=object$par[1:nit]; vi=object$par[(nit+1):(2*nit)]; ter=object$par[(2*nit+1):(3*nit)]
  oA=object$par[3*nit+1]; oV=object$par[3*nit+2]
  model=if(object$model==1) "D" else "Q"
  miss=is.na(object$rt) | is.na(object$score)
  out=lapply(1:nsim,function(k){
    gamma=exp(oA*rnorm(N))
    theta=if(model=="D") oV*rnorm(N) else exp(oV*rnorm(N))
    s=simdiff(N,nit,ai=ai,vi=vi,gamma=gamma,theta=theta,ter=ter,model=model,...)
    if(keep.missing){ s$rt[miss]=NA; s$x[miss]=NA }
    s
  })
  attr(out,"seed")=seed
  out
}
