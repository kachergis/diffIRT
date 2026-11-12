# Item fit of the joint distribution of responses and response times.
#
# For each item, the responses of the persons who answered it are classified into 2 x bins
# cells: response (0/1) crossed with RT bins whose edges are the deciles (for bins = 10) of
# the item's predicted marginal RT distribution. The observed cell proportions p are compared
# with the model's pi(theta-hat), integrated over the person distribution.
#
# Covariance of the residuals p - pi(theta-hat). For persons answering item j,
# Cov(p_j) = (diag(pi_j) - pi_j pi_j') / n_j; for two items the cells overlap in the n_jk
# persons who answered both, Cov(p_j, p_k) = n_jk / (n_j n_k) (pi_jk - pi_j pi_k'), with pi_jk
# the joint cell probabilities. Because theta-hat is the maximum likelihood estimate,
# Cov(p, theta-hat) = Delta V, with Delta = d pi / d theta and V = vcov(object), so
#   Cov(p - pi(theta-hat)) = Cov(p) - Delta V Delta'
# (the same argument as for limited-information tests, here with the full-information MLE).
# The statistic r' Cov^- r is asymptotically chi-square with df = the rank of the covariance
# (2 x bins - 1 per item). Cells defined by estimated quantiles have the same limiting
# distribution as cells fixed at the true quantiles (Moore & Spruill, 1975). Missing
# responses are assumed missing completely at random.

RTFit=function(object,bins=10,control=list()){
  if(object$conv!=0) stop("The solution in 'object' has not converged.")
  if(bins<2) stop("'bins' should be at least 2.")
  con=list(nq=21,panels=12,h=1e-4); con[names(control)]=control
  nit=object$nit; J=nit; N=object$N; K=bins
  par=object$par.log; constrain=object$constrain; npars=object$npars
  x=object$score; rt=object$rt
  obs=!is.na(x) & !is.na(rt)

  V=fitVcov(object)

  # bin edges: predicted marginal quantiles of each item's RT
  edges=lapply(1:J,function(j){
    mx=2*max(rt[,j],na.rm=TRUE)
    qdiffIRT((1:(K-1))/K,par[j],par[J+j],par[2*J+j],par[3*J+1],par[3*J+2],mx,
             object$A,object$W,object$model,object$control$eps)
  })

  gq=statmod::gauss.quad(con$nq,kind="hermite"); A=gq$nodes*sqrt(2); W=gq$weights/sqrt(pi)
  w=as.vector(outer(W,W))
  gl=statmod::gauss.quad(8)
  cells=function(p,j) rtCellNodes(p,j,J,edges[[j]],object$model,object$control$eps,A,gl,con$panels)

  M=lapply(1:J,function(j) cells(par,j))                  # nodes x 2K, per item
  pis=lapply(M,function(m) colSums(w*m))

  # observed proportions
  cnt=lapply(1:J,function(j){
    o=obs[,j]; b=findInterval(rt[o,j],edges[[j]])+1
    tabulate(x[o,j]*K+b,2*K)
  })
  n=colSums(obs)
  pobs=lapply(1:J,function(j) cnt[[j]]/n[j])

  # derivatives of the cell probabilities with respect to the free parameters
  h=con$h
  D=matrix(0,J*2*K,npars)
  for(q in 1:npars){
    k=which(constrain==q)
    its=if(any(k>3*J)) 1:J else unique(((k-1)%%J)+1)
    pp=par; pm=par; pp[k]=pp[k]+h; pm[k]=pm[k]-h
    for(j in its){
      rows=(j-1)*2*K+1:(2*K)
      D[rows,q]=(colSums(w*cells(pp,j))-colSums(w*cells(pm,j)))/(2*h)
    }
  }

  keep=as.vector(sapply(1:J,function(j) (j-1)*2*K+1:(2*K-1)))  # drop each item's last cell
  r=unlist(pobs)-unlist(pis)
  nn=crossprod(obs+0)
  S=matrix(0,J*2*K,J*2*K)
  for(j in 1:J) for(k in j:J){
    rj=(j-1)*2*K+1:(2*K); rk=(k-1)*2*K+1:(2*K)
    cv=if(j==k) diag(pis[[j]])-tcrossprod(pis[[j]]) else crossprod(w*M[[j]],M[[k]])-tcrossprod(pis[[j]],pis[[k]])
    S[rj,rk]=cv*nn[j,k]/(n[j]*n[k])
    S[rk,rj]=t(S[rj,rk])
  }
  S=S-D%*%V%*%t(D)

  wald=function(idx){
    e=eigen(S[idx,idx],symmetric=TRUE)
    ok=e$values>1e-8*max(e$values)
    u=crossprod(e$vectors[,ok,drop=FALSE],r[idx])
    c(stat=sum(u^2/e$values[ok]),df=sum(ok))
  }
  it=t(sapply(1:J,function(j) wald(intersect(keep,(j-1)*2*K+1:(2*K)))))
  items=data.frame(n=n,stat=it[,"stat"],df=it[,"df"],p=pchisq(it[,"stat"],it[,"df"],lower.tail=FALSE),
                   row.names=paste("item",1:J))
  ov=wald(keep)

  lab=c(paste0("x=0 bin",1:K),paste0("x=1 bin",1:K))
  tab=lapply(1:J,function(j) data.frame(obs=pobs[[j]],exp=pis[[j]],row.names=lab))
  names(tab)=rownames(items)
  res=list(items=items,overall=c(stat=unname(ov["stat"]),df=unname(ov["df"]),
                                 p=pchisq(ov["stat"],ov["df"],lower.tail=FALSE)[[1]]),
           cells=tab,edges=edges,bins=K,model=object$model)
  class(res)="RTFit"
  res
}

# Cell probabilities of item j at each node of a Q x Q product grid over the prior
# (rows: node a + Q (b - 1); columns: x = 0 bins 1..K, then x = 1 bins 1..K). The bins
# below the last edge are integrated with 8-point Gauss-Legendre on 'panels' panels per
# bin (quadratically spaced in the first bin, where the density rises from 0); the last
# bin is P(x | node) minus the others.
rtCellNodes=function(par,j,J,e,model,eps,A,gl,panels){
  Q=length(A); K=length(e)+1
  ai=par[j]; vi=par[J+j]; ter=exp(par[2*J+j])
  zA=sqrt(exp(par[3*J+1]))*A; zV=sqrt(exp(par[3*J+2]))*A
  la=zA-ai; a=exp(la)
  v=if(model==1) zV-vi else exp(zV-vi)
  ed=pmax(c(0,e-ter),0)                                   # decision-time bin edges
  s=seq(0,1,length=panels+1)
  pts=lapply(1:(K-1),function(k){
    pe=ed[k]+(ed[k+1]-ed[k])*(if(k==1) s^2 else s)
    hw=diff(pe)/2
    list(t=as.vector(outer(gl$nodes+1,hw)+rep(pe[-length(pe)],each=8)),
         w=as.vector(outer(gl$weights,hw)))
  })
  t=unlist(lapply(pts,`[[`,"t")); gw=unlist(lapply(pts,`[[`,"w"))
  bin=rep(1:(K-1),each=8*panels)
  B=matrix(.C("ddiffBase",as.double(t),as.integer(length(t)),as.double(la),as.integer(Q),
              as.double(eps),out=double(length(t)*Q),PACKAGE="diffIRT")$out,length(t),Q)
  out=matrix(0,Q*Q,2*K)
  vt=outer(t,v^2)/2                                       # T x Q (drift nodes)
  for(qa in 1:Q){
    rows=qa+Q*(0:(Q-1))                                   # nodes (qa, b) for all b
    p1=1/(1+exp(-a[qa]*v))
    for(c in 0:1){
      f=exp(B[,qa]+outer(rep(1,length(t)),v*a[qa]*(c-.5))-vt)
      f[!is.finite(f)]=0
      bi=rowsum(gw*f,bin,reorder=TRUE)                    # (K-1) x Q
      pc=if(c==1) p1 else 1-p1
      out[rows,c*K+1:(K-1)]=t(bi)
      out[rows,c*K+K]=pc-colSums(bi)
    }
  }
  out
}

# Covariance matrix of the free-parameter estimates; computed from the fit if it was
# not fitted with se = TRUE
fitVcov=function(object){
  if(!is.null(object$hessian)) return(vcov(object))
  rt=object$rt; x=object$score
  rt[is.na(rt)]=-999; rt[is.na(x)]=-999; x[is.na(x)]=-999
  pars=object$par.log[match(1:object$npars,object$constrain)]
  con=object$control
  H=optimHess(pars,likelihood,deriv,rt=rt,x=x,N=object$N,nit=object$nit,model=object$model,A=object$A,
              W=object$W,nq=length(object$A),eps=con$eps,delta=con$delta,start=object$par.log,
              constrain=object$constrain,grad=if(is.null(con$grad)) "analytic" else con$grad,aq=object$aq)
  solve(.5*H)
}

print.RTFit=function(x,digits=3,...){
  cat("Item fit of responses and response times (",x$bins," RT bins x response)\n\n",sep="")
  it=x$items; it$stat=round(it$stat,2); it$p=round(it$p,digits)
  print(it)
  cat("\nOverall: stat = ",round(x$overall["stat"],2),", df = ",x$overall["df"],", p = ",
      format.pval(x$overall["p"],digits=digits),"\n",sep="")
  invisible(x)
}
