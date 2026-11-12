RespFit=function(object,order=2){
  # Limited-information goodness-of-fit test Mr (Maydeu-Olivares & Joe, 2005) for the responses.
  #
  # Mr only needs the marginal moments up to 'order' (pi_S = P(all items in S correct),
  # |S| <= order), their asymptotic covariance and their derivatives. With binary items
  # x_j^2 = x_j, so Cov(x_S, x_T) = pi_{S u T} - pi_S pi_T: every quantity needs joint
  # probabilities of at most 2*order items, computed directly by quadrature. (Earlier
  # versions enumerated all 2^nit response patterns, which was infeasible beyond ~14 items.)
  if(object$conv!=0) stop("The solution in 'object' has not converged.")
  if(object$nr_fail!=0) cat("Warning: The solution provided in 'object' might be instable.\n")
  if(order==1) stop("The Mr test is not possible with only first-order moments, please use at least order=2")

  nit=object$nit
  N=object$N
  model=object$model
  W=object$W
  A=object$A
  x=object$score
  constrain=object$constrain
  if(order>nit) stop("'order' cannot exceed the number of items.")

  # the moments: item sets S (one per row of fitItems, without its leading all-zero row),
  # in the order: all single items, then all pairs, ...
  X=fitItems(nit,order)[-1,,drop=FALSE]
  sets=lapply(1:nrow(X),function(r) which(X[r,]==1))
  nm=length(sets)

  # Observed moments from pairwise-available data: p_S is the proportion of persons with all
  # items in S correct, among the n_S persons who answered all items in S. (Versions <= 1.5
  # divided by N, so a missing response counted as incorrect.) Under MCAR the persons behind
  # p_S and p_T overlap in the n_{S u T} who answered all of S u T, so
  # Cov(p_S, p_T) = n_{S u T} / (n_S n_T) * (pi_{S u T} - pi_S pi_T); with no missing data this
  # is the usual Xi / N.
  obs=!is.na(x)
  x0=x; x0[!obs]=0
  R=sapply(sets,function(S) rowSums(obs[,S,drop=FALSE])==length(S))          # N x moments: answered all of S
  C=sapply(sets,function(S) rowSums(x0[,S,drop=FALSE])==length(S))           # N x moments: all of S correct
  R=matrix(R,N); C=matrix(C,N)
  n=colSums(R)
  if(any(n<2)) stop("Some items (or pairs of items) have fewer than 2 joint responses; Mr cannot be computed.")
  o=colSums(C)/n
  nn=crossprod(R+0)/outer(n,n)                 # n_{S u T} / (n_S n_T)

  # probability correct for each item at each node of the (boundary x drift) product grid
  nq=length(A)
  Wq=as.vector(outer(W,W))                     # node (n1, n2) -> index n1 + nq*(n2-1)
  probs=function(th){
    ai=th[1:nit]; vi=th[(nit+1):(2*nit)]; sA=sqrt(exp(th[2*nit+1])); sV=sqrt(exp(th[2*nit+2]))
    gam=rep(exp(A*sA),nq)                      # person boundary at node n1
    lat=rep(A*sV,each=nq)                      # person drift latent at node n2
    if(model==1) logit=outer(gam,exp(ai),"/")*outer(lat,vi,"-")
    else logit=outer(gam,exp(ai),"/")*outer(exp(lat),exp(vi),"/")
    logit[logit<(-709)]=-709
    1/(1+exp(-logit))                          # nodes x items
  }
  moments=function(p) sapply(sets,function(S) sum(Wq*apply(p[,S,drop=FALSE],1,prod)))

  th=object$par.log[c(1:(2*nit),3*nit+1,3*nit+2)]
  p=probs(th)
  M=sapply(sets,function(S) apply(p[,S,drop=FALSE],1,prod))   # nodes x moments
  pi.vec=colSums(Wq*M)

  # E[x_S x_T] = pi_{S u T}: products of the node-wise moments are right unless S and T share items
  E2=crossprod(M*Wq,M)
  for(a in 1:nm) for(b in a:nm){
    sh=intersect(sets[[a]],sets[[b]])
    if(length(sh)){
      E2[a,b]=E2[b,a]=sum(Wq*apply(p[,union(sets[[a]],sets[[b]]),drop=FALSE],1,prod))
    }
  }
  k=(E2-pi.vec%*%t(pi.vec))*nn                  # covariance of the observed moments (Xi_r / N without missing data)

  # Derivatives of the moments with respect to the free parameters (mapped through
  # 'constrain'). Responses alone do not identify every parameter: in the D model the
  # scale of the drift is free (multiplying theta, v, a and omega_theta by c leaves all
  # probabilities unchanged); in the Q model accuracy depends on a[i]*v[i] and on
  # omega_gamma^2 + omega_theta^2 only. Mr is therefore computed with an orthonormal basis
  # of the column space of the derivative matrix, and df = #moments - its rank.
  # (In versions <= 1.5 the Q-model derivatives were wrong: fitDeriv tested
  # toupper(model)=="Q", which is never true for the numeric model code 2, so v* stayed on
  # the log scale. In the D model the derivatives were right, but the free drift scale was
  # treated as identified, and the df subtracted 2 more than the number of parameters.)
  h=1e-5
  dfull=sapply(seq_along(th),function(j){ tp=th; tm=th; tp[j]=tp[j]+h; tm[j]=tm[j]-h
    (moments(probs(tp))-moments(probs(tm)))/(2*h) })
  cpos=constrain[c(1:(2*nit),3*nit+1,3*nit+2)]  # free-parameter number of each column (0 = fixed)
  free=sort(unique(cpos[cpos!=0]))
  d=sapply(free,function(q) rowSums(dfull[,cpos==q,drop=FALSE]))
  d=matrix(d,nm)
  d=d[,colSums(d^2)>0,drop=FALSE]
  sv=svd(sweep(d,2,sqrt(colSums(d^2)),"/"))
  rk=sum(sv$d>1e-5*sv$d[1])            # true directions are >1e-2; quadrature artifacts ~1e-7 or less
  U=sv$u[,1:rk,drop=FALSE]

  Mr=-999
  df=-999
  Cr=matrix(NA,nm,nm)
  ki=try(solve(k),silent=TRUE)
  if(inherits(ki,"try-error")){
    cat("error: the covariance matrix of the moments is singular\n")
  } else {
    Cr=ki-ki%*%U%*%solve(t(U)%*%ki%*%U)%*%t(U)%*%ki
    res=o-pi.vec
    Mr=as.numeric(t(res)%*%Cr%*%res)
    df=nm-rk
  }

  nms=sapply(1:nrow(X),function(r) paste(X[r,],collapse=""))
  nrms=(o-pi.vec)*sqrt(diag(Cr))
  mat=data.frame(pred=round(pi.vec,3),obs=round(o,3),Z=round(nrms,3),row.names=nms)
  respfitRES=list(model=model,Z=mat,Mr=Mr,df=df,order=order,nit=nit)
  class(respfitRES)<- "RespFit"
  return(respfitRES)
}
