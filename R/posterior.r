# Posterior quantities of the person parameters for factest(type = "EAP" / "PV").
#
# The latent variables are zA = log gamma_p ~ N(0, omega_gamma^2) and zV ~ N(0, omega_theta^2),
# with theta_p = zV (D model) or exp(zV) (Q model). LLpoints (src/diffIRT.c) gives each
# person's log-likelihood at any (zA, zV); the prior is added here.

# log-likelihood at points: zA, zV are N x M matrices (row = person); returns N x M
llPoints=function(par,rt,x,N,nit,model,eps,zA,zV){
  M=ncol(zA)
  r=.C("LLpoints",as.double(rt),as.double(x),as.integer(N),as.integer(nit),as.double(par[1:nit]),
       as.double(par[(nit+1):(2*nit)]),as.double(par[(2*nit+1):(3*nit)]),as.integer(model),as.double(eps),
       as.integer(M),as.double(t(zA)),as.double(t(zV)),ou=double(N*M),PACKAGE="diffIRT")
  matrix(r$ou,N,M,byrow=TRUE)
}

# log prior density of z under N(0, s2), up to a constant; 0 when the variance is 0
logPrior=function(z,s2) if(s2>0) -z^2/(2*s2) else 0*z

# per-person posterior centre and spread (log gamma, drift latent). For the fitted parameters
# these are the fit's adaptive grids (object$aq). Otherwise (other parameters, e.g. drawn for
# plausible values, or a non-adaptive fit) they are recomputed with two passes of aqMoments,
# starting from the fit's adaptive grids when there are any.
postCentres=function(object,par,rt,x,fitted=TRUE){
  nit=object$nit; P=3*nit+2
  if(fitted && !is.null(object$aq)) return(object$aq)
  mom=function(aq) aqMoments(par,rt,x,object$N,nit,object$model,object$A,object$W,length(object$A),
                             object$control$eps,par,1:P,aq=aq)
  mom(mom(object$aq))
}

# weights of each person's posterior on an adaptive Gauss-Hermite product grid with Q nodes
# per dimension, centred and scaled on the person's posterior (aq)
postGrid=function(object,par,rt,x,aq,Q){
  N=object$N; nit=object$nit
  gq=statmod::gauss.quad(Q,kind="hermite"); A=gq$nodes*sqrt(2); lw=log(gq$weights)+A^2/2
  s2A=exp(par[3*nit+1]); s2V=exp(par[3*nit+2])
  ia=rep(1:Q,Q); ib=rep(1:Q,each=Q)                 # product grid, Q^2 points
  zA=if(s2A>0) aq$muA+outer(aq$sA,A[ia]) else matrix(0,N,Q*Q)
  zV=if(s2V>0) aq$muV+outer(aq$sV,A[ib]) else matrix(0,N,Q*Q)
  lw=llPoints(par,rt,x,N,nit,object$model,object$control$eps,zA,zV)+logPrior(zA,s2A)+logPrior(zV,s2V)+
     rep(lw[ia]+lw[ib],each=N)
  m=apply(lw,1,max)
  w=exp(lw-m); w=w/rowSums(w)
  list(zA=zA,zV=zV,w=w,ok=is.finite(m))
}

# Expected a posteriori scores and posterior SDs on the reported scale
eapScores=function(object,par,rt,x,aq,Q){
  g=postGrid(object,par,rt,x,aq,Q)
  tV=if(object$model==2) exp else identity
  gam=exp(g$zA); th=tV(g$zV)
  mA=rowSums(g$w*gam); mV=rowSums(g$w*th)
  out=cbind(mA,sqrt(pmax(rowSums(g$w*gam^2)-mA^2,0)),mV,sqrt(pmax(rowSums(g$w*th^2)-mV^2,0)))
  out[!g$ok,]=NA
  out
}

# One plausible value per person: importance resampling from K draws of a bivariate t(5)
# proposal centred on the person's posterior with 1.5 times its SD
pvDraw=function(object,par,rt,x,aq,K){
  N=object$N; nit=object$nit
  s2A=exp(par[3*nit+1]); s2V=exp(par[3*nit+2])
  df=5; inflate=1.5
  tA=matrix(stats::rt(N*K,df),N); tV=matrix(stats::rt(N*K,df),N)   # (rt is the data here)
  zA=if(s2A>0) aq$muA+inflate*aq$sA*tA else matrix(0,N,K)
  zV=if(s2V>0) aq$muV+inflate*aq$sV*tV else matrix(0,N,K)
  lq=if(s2A>0) stats::dt(tA,df,log=TRUE) else 0
  lq=lq+if(s2V>0) stats::dt(tV,df,log=TRUE) else 0
  lw=llPoints(par,rt,x,N,nit,object$model,object$control$eps,zA,zV)+logPrior(zA,s2A)+logPrior(zV,s2V)-lq
  m=apply(lw,1,max)
  w=exp(lw-m)
  ess=rowSums(w)^2/rowSums(w^2)
  k=apply(w,1,function(p) if(all(is.finite(p)) && sum(p)>0) sample.int(K,1,prob=p) else NA)
  i=cbind(1:N,k)
  list(zA=zA[i],zV=zV[i],ess=ess)
}
