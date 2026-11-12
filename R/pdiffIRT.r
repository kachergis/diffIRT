# Predicted (marginal) RT distribution of one item, both responses combined.
#
# cdfGrid() integrates the density mdiffIRT() once over [Ter, maxRT] on 'panels' panels
# with 8-point Gauss-Legendre in each panel. The panel edges are Ter + (maxRT-Ter)*s^2
# for equally spaced s, so the panels are narrow where the density rises steeply just
# after Ter and wide in the smooth right tail. This gives the CDF at the panel edges
# (and the density there). Between edges the CDF is the cubic Hermite interpolant of
# these values and slopes, accurate to O(h^4). pdiffIRT() and qdiffIRT() evaluate and
# invert this interpolant, so all quantiles of an item cost one pass over the grid.
# (Versions <= 1.5 called integrate() for every RT inside uniroot() for every quantile.)
# If the CDF at maxRT is below 'pmax', the range is doubled until it is not.
cdfGrid=function(ai,vi,ter,sd2A,sd2V,maxRT,A,W,model,eps,panels=500,pmax=.999){
  gl=statmod::gauss.quad(8)                       # nodes/weights on [-1,1]
  lo=exp(ter)
  if(!(maxRT>lo)) stop("'maxRT' should exceed the non-decision time.")
  repeat{
    e=lo+(maxRT-lo)*seq(0,1,length=panels+1)^2; h=diff(e)
    t=as.vector(outer((gl$nodes+1)/2,h)+rep(e[-(panels+1)],each=8))
    d=mdiffIRT(c(t,e),ai,vi,ter,sd2A,sd2V,A,W,model,eps)
    d[is.na(d)]=0
    inc=colSums(matrix(d[seq_along(t)],8)*gl$weights)*h/2
    Fe=c(0,cumsum(inc))
    fe=d[length(t)+seq_along(e)]
    if(Fe[panels+1]>=pmax || maxRT>1e4*lo) break
    maxRT=lo+2*(maxRT-lo)
  }
  list(e=e,F=Fe,f=fe,h=h)
}

# cubic Hermite CDF at y (vector), given the grid
cdfEval=function(g,y){
  n=length(g$e)
  k=pmin(pmax(findInterval(y,g$e),1),n-1)
  hk=g$h[k]
  s=pmin(pmax((y-g$e[k])/hk,0),1)
  h00=2*s^3-3*s^2+1; h10=s^3-2*s^2+s; h01=-2*s^3+3*s^2; h11=s^3-s^2
  out=h00*g$F[k]+h10*hk*g$f[k]+h01*g$F[k+1]+h11*hk*g$f[k+1]
  out[y<=g$e[1]]=0
  out[y>=g$e[n]]=g$F[n]
  out
}

pdiffIRT=function(y,ai,vi,ter,sd2A,sd2V,A,W,model,eps,maxRT=max(y)){
  g=cdfGrid(ai,vi,ter,sd2A,sd2V,max(maxRT,exp(ter)+1e-3),A,W,model,eps,pmax=0)
  cdfEval(g,y)
}
