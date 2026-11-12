# Quantiles of the predicted RT distribution of one item (see cdfGrid in pdiffIRT.r):
# locate each probability's panel, then bisect the cubic Hermite CDF within it
# (vectorized over p; 50 halvings of one panel).
qdiffIRT=function(p,ai,vi,ter,sd2A,sd2V,maxRT,A,W,model,eps,panels=500){
  g=cdfGrid(ai,vi,ter,sd2A,sd2V,maxRT,A,W,model,eps,panels=panels,pmax=max(p))
  if(max(p)>g$F[length(g$F)]) stop("The predicted CDF does not reach the requested probabilities; increase 'maxRT'.")
  k=pmin(findInterval(p,g$F,left.open=TRUE),length(g$e)-1)
  lo=g$e[k]; hi=g$e[k+1]
  for(it in 1:50){
    mid=(lo+hi)/2
    below=cdfEval(g,mid)<p
    lo[below]=mid[below]; hi[!below]=mid[!below]
  }
  (lo+hi)/2
}
