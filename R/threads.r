# Set the number of OpenMP threads for later parallel regions and return the
# previous setting. NULL or NA leaves the OpenMP default (OMP_NUM_THREADS, or all
# cores) unchanged. Callers restore the previous value with on.exit().
setThreads=function(n){
  if(is.null(n) || (length(n)==1 && is.na(n))) n=0L
  if(length(n)!=1 || !is.numeric(n) || n<0 || n!=round(n))
    stop("control$threads should be a single positive integer (or NULL for the OpenMP default)")
  .C("diffIRT_threads",as.integer(n),integer(1),PACKAGE="diffIRT")[[2]]
}
