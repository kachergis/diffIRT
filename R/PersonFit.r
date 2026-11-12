# Person fit: each person's marginal log-likelihood l_i = log p(x_i, rt_i), compared with
# its distribution under the fitted model for a person with the same items answered.
#
# The reference distribution is obtained by Monte Carlo: nsim persons are simulated from the
# fitted model once (person parameters drawn from the population distribution), and for each
# pattern of answered items their log-likelihood is computed over those items only. Because
# the observed and simulated statistics are the same function of the data (log-likelihood on
# a fixed nq x nq Gauss-Hermite grid), the p-values are calibrated whatever the accuracy of
# that grid, apart from the estimation of the item parameters (small for large N).
# Aberrant persons (e.g. rapid guessing, or responses inconsistent with the RTs) have a low
# l_i: p is the lower-tail Monte Carlo p-value, lz = (l_i - mean) / sd of the reference.
# Item parameters estimated with the aberrant persons included can mask them: rapid guesses
# pull Ter down, after which they look like ordinary fast responses. (A purification step,
# refitting without the flagged persons, gained little power and raised the type I error,
# so it is not offered.)

PersonFit=function(object,nsim=1000,seed=NULL,control=list()){
  if(object$conv!=0) stop("The solution in 'object' has not converged.")
  con=list(nq=11,threads=object$control$threads); con[names(control)]=control
  if(!is.null(con$threads)){ oldThreads=setThreads(con$threads); on.exit(setThreads(oldThreads),add=TRUE) }
  if(!is.null(seed)) set.seed(seed)
  out=list(persons=personStats(object,object$par.log,nsim,con$nq),nsim=nsim,nq=con$nq)
  class(out)="PersonFit"
  out
}

# statistics of all persons in 'object' under parameters 'par' (internal scale)
personStats=function(object,par,nsim,nq){
  J=object$nit; N=object$N
  gq=statmod::gauss.quad(nq,kind="hermite"); A=gq$nodes*sqrt(2); W=gq$weights/sqrt(pi)
  ll=function(rt,x){
    rt[is.na(rt)]=-999; rt[is.na(x)]=-999; x[is.na(x)]=-999
    r=.C("LLdiff",as.double(rt),as.double(x),as.integer(nrow(rt)),as.integer(J),as.double(par[1:J]),
         as.double(par[J+1:J]),as.double(par[3*J+1]),as.double(par[3*J+2]),as.double(par[2*J+1:J]),
         as.integer(object$model),as.double(A),as.double(W),as.integer(nq),as.double(object$control$eps),
         as.double(-1),LL=double(nrow(rt)),PACKAGE="diffIRT")
    log(r$LL)
  }
  l=ll(object$rt,object$score)
  sim=simPersons(par,J,object$model,nsim)        # simulated persons, all items answered
  obs=!is.na(object$rt) & !is.na(object$score)
  key=apply(obs,1,function(o) paste(which(o),collapse=","))
  mu=sdv=p=numeric(N)
  for(k in unique(key)){
    i=which(key==k); o=obs[i[1],]
    rt=sim$rt; x=sim$x; rt[,!o]=NA; x[,!o]=NA
    ls=ll(rt,x)
    mu[i]=mean(ls); sdv[i]=stats::sd(ls)
    p[i]=sapply(l[i],function(li) (1+sum(ls<=li))/(nsim+1))
  }
  data.frame(ll=l,lz=(l-mu)/sdv,p=p,items=rowSums(obs))
}

# nsim persons simulated with all items answered, under parameters 'par' (internal scale)
simPersons=function(par,J,model,nsim){
  sA=sqrt(exp(par[3*J+1])); sV=sqrt(exp(par[3*J+2]))
  gamma=exp(sA*stats::rnorm(nsim))
  theta=if(model==1) sV*stats::rnorm(nsim) else exp(sV*stats::rnorm(nsim))
  simdiff(nsim,J,ai=exp(par[1:J]),vi=if(model==1) par[J+1:J] else exp(par[J+1:J]),gamma=gamma,theta=theta,
          ter=exp(par[2*J+1:J]),model=if(model==1) "D" else "Q")
}

print.PersonFit=function(x,alpha=.05,...){
  pf=x$persons
  cat("Person fit: marginal log-likelihood against ",x$nsim," simulated persons per pattern of answered items\n\n",sep="")
  cat(sum(pf$p<alpha)," of ",nrow(pf)," persons (",round(100*mean(pf$p<alpha),1),"%) have p < ",alpha,
      " (expected under the model: ",100*alpha,"%)\n",sep="")
  worst=order(pf$p,pf$lz)[1:min(10,nrow(pf))]
  cat("\nLowest p-values:\n")
  print(round(pf[worst,],3))
  invisible(x)
}
