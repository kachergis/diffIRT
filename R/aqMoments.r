# Posterior mean and SD of each subject's log boundary (gamma) and drift latent,
# computed on the current quadrature grid (aq = NULL: the non-adaptive grid).
# Used to re-centre the adaptive Gauss-Hermite grid between optimizer runs.
aqMoments=function(pars,rt,x,N,nit,model,A,W,nq,eps,start,constrain,aq=NULL){
    par=rep(999,3*nit+2)
    par[constrain!=0]=pars[constrain]
    par[constrain==0]=start[constrain==0]
    sd2A=par[3*nit+1]
    sd2V=par[3*nit+2]
    if(is.null(aq)) aq=list(muA=rep(0,N),sA=rep(sqrt(exp(sd2A)),N),muV=rep(0,N),sV=rep(sqrt(exp(sd2V)),N))
    r=.C("AQmoments",as.double(rt),as.double(x),as.integer(N),as.integer(nit),as.double(par[1:nit]),
         as.double(par[(nit+1):(2*nit)]),as.double(sd2A),as.double(sd2V),as.double(par[(2*nit+1):(3*nit)]),as.integer(model),
         as.double(A),as.double(W),as.integer(nq),as.double(eps),
         as.double(aq$muA),as.double(aq$sA),as.double(aq$muV),as.double(aq$sV),as.double(-1),double(4*N))
    m=matrix(r[[20]],4)
    # keep the grid from collapsing (one node holding all the weight)
    list(muA=m[1,],sA=pmax(m[2,],1e-3),muV=m[3,],sV=pmax(m[4,],1e-3),LL=r[[19]])
}
