#!/usr/bin/env Rscript
# benchmark.R
# Compares original vs optimised LLdiff by calling the compiled C directly.
# Run from the repo root: Rscript benchmark.R

suppressPackageStartupMessages({
  library(statmod)   # gauss.quad
  library(microbenchmark)
})

# --------------------------------------------------------------------------
# 1. Compile both versions to shared libraries
# --------------------------------------------------------------------------
src_orig <- "src_orig/diffIRT.c"   # original (from CRAN / your repo)
src_opt  <- "src/diffIRT.c"        # optimised version

compile_so <- function(src, out, openmp = FALSE) {
  # NB: "-fopenmp" must reach the *compiler*. Passing it as an extra argument to
  # R CMD SHLIB only adds it to the link line, so _OPENMP was never defined and
  # the earlier "optimised" build silently ran serially. PKG_CFLAGS fixes that.
  omp_env <- if (openmp) "PKG_CFLAGS=-fopenmp PKG_LIBS=-fopenmp" else ""
  cmd <- sprintf(
    "%s R CMD SHLIB -o %s %s 2>&1",
    omp_env, out, src
  )
  res <- system(cmd, intern = TRUE)
  if (!file.exists(out)) {
    cat("Compilation failed:\n"); cat(res, sep="\n"); stop()
  }
  invisible(res)
}

dir.create("src_orig", showWarnings = FALSE)
# Write the original C if not already present (paste from document or copy)
if (!file.exists("src_orig/diffIRT.c"))
  stop("Please place the ORIGINAL diffIRT.c in src_orig/diffIRT.c")

cat("Compiling original (no OpenMP)...\n")
compile_so("src_orig/diffIRT.c", "src_orig/diffIRT_orig.so", openmp = FALSE)

cat("Compiling optimised (with OpenMP if available)...\n")
compile_so("src/diffIRT.c", "src/diffIRT_opt.so", openmp = TRUE)

dyn.load("src_orig/diffIRT_orig.so")
dyn.load("src/diffIRT_opt.so")

# --------------------------------------------------------------------------
# 2. Simulate data
# --------------------------------------------------------------------------
set.seed(42)
N   <- 500
nit <- 20

# True parameters
ai  <- exp(runif(nit, -1, 0))
vi  <- runif(nit, -1.5, 0.5)
ter <- runif(nit, 0.3, 0.8)
sd2A <- log(0.3^2)   # log-scale as used internally
sd2V <- log(0.5^2)

# Quick simulation (rejection sampler from simdiff would be ideal;
# use Wald approximation for benchmark purposes)
gamma_p <- rlnorm(N, 0, exp(sd2A/2))
theta_p <- rnorm(N, 0, exp(sd2V/2))

rt <- matrix(NA, N, nit)
x  <- matrix(NA, N, nit)
for (j in 1:nit) {
  for (p in 1:N) {
    a  <- gamma_p[p] / ai[j]
    mu <- theta_p[p] - vi[j]
    # Wald mean + noise as proxy RT
    rt[p, j] <- abs(a / (mu + 1e-6)) + ter[j] + abs(rnorm(1, 0, 0.1))
    x[p, j]  <- as.integer(runif(1) < 1 / (1 + exp(-mu * a)))
  }
}
rt[rt < 0.05] <- 0.05

# --------------------------------------------------------------------------
# 3. Set up quadrature
# --------------------------------------------------------------------------
nq <- 7
Quad <- gauss.quad(nq, kind = "hermite")
W <- Quad$weights / sqrt(pi)
A <- Quad$nodes  * sqrt(2)
eps <- 0.01

# Internal log-scale parameters (as passed by diffIRT's R wrapper)
ai_log  <- log(ai)
ter_log <- log(ter)

# Use log-scale sd2 (already in log scale above)
# vi is not log-transformed for D model

# Sentinel for missing: -999
rt_c <- rt;  rt_c[is.na(rt_c)] <- -999
x_c  <- x;   x_c[is.na(x_c)]  <- -999
model <- 1L   # D-diffusion

# Flatten to vectors as .C expects column-major
rt_vec <- as.double(rt_c)
x_vec  <- as.double(x_c)
LL_vec <- as.double(rep(0, N))

call_orig <- function() {
  .C("LLdiff",
     rt_vec, x_vec,
     as.integer(N), as.integer(nit),
     as.double(ai_log), as.double(vi),
     as.double(sd2A), as.double(sd2V),
     as.double(ter_log),
     as.integer(model),
     as.double(A), as.double(W), as.integer(nq),
     as.double(eps),
     as.double(-1), as.double(rep(0, N)),
     PACKAGE = "diffIRT_orig")
}

call_opt <- function() {
  .C("LLdiff",
     rt_vec, x_vec,
     as.integer(N), as.integer(nit),
     as.double(ai_log), as.double(vi),
     as.double(sd2A), as.double(sd2V),
     as.double(ter_log),
     as.integer(model),
     as.double(A), as.double(W), as.integer(nq),
     as.double(eps),
     as.double(-1), as.double(rep(0, N)),
     PACKAGE = "diffIRT_opt")
}

# --------------------------------------------------------------------------
# 4. Sanity check: results should be numerically identical
# --------------------------------------------------------------------------
res_orig <- call_orig()
res_opt  <- call_opt()

cat(sprintf("\nLL original : %.6f\n", res_orig[[15]]))
cat(sprintf("LL optimised: %.6f\n", res_opt[[15]]))
cat(sprintf("Max abs diff in subject LLs: %.2e\n\n",
            max(abs(res_orig[[16]] - res_opt[[16]]))))

# --------------------------------------------------------------------------
# 5. Benchmark
# --------------------------------------------------------------------------
cat("Benchmarking (100 reps)...\n")
mb <- microbenchmark(
  original  = call_orig(),
  optimised = call_opt(),
  times = 100L
)
print(mb)
cat("\nSpeedup (median): ",
    round(median(mb$time[mb$expr=="original"]) /
          median(mb$time[mb$expr=="optimised"]), 2), "x\n")

# --------------------------------------------------------------------------
# 6. Scaling test: vary N
# --------------------------------------------------------------------------
cat("\n--- Scaling with N (nit=20, nq=7) ---\n")

# Generate larger dataset for scaling tests
set.seed(123)
N_max <- 2000
gamma_p_large <- rlnorm(N_max, 0, exp(sd2A/2))
theta_p_large <- rnorm(N_max, 0, exp(sd2V/2))

rt_large <- matrix(NA, N_max, nit)
x_large  <- matrix(NA, N_max, nit)
for (j in 1:nit) {
  for (p in 1:N_max) {
    a  <- gamma_p_large[p] / ai[j]
    mu <- theta_p_large[p] - vi[j]
    rt_large[p, j] <- abs(a / (mu + 1e-6)) + ter[j] + abs(rnorm(1, 0, 0.1))
    x_large[p, j]  <- as.integer(runif(1) < 1 / (1 + exp(-mu * a)))
  }
}
rt_large[rt_large < 0.05] <- 0.05

rt_c_large <- rt_large; rt_c_large[is.na(rt_c_large)] <- -999
x_c_large  <- x_large;  x_c_large[is.na(x_c_large)]   <- -999

for (Ns in c(100, 250, 500, 1000, 2000)) {
  rt_s <- as.double(rt_c_large[1:Ns, ])
  x_s  <- as.double(x_c_large[1:Ns, ])
  Ni   <- as.integer(Ns)
  t_o <- median(microbenchmark(
    .C("LLdiff", rt_s, x_s, Ni, as.integer(nit),
       as.double(ai_log), as.double(vi),
       as.double(sd2A), as.double(sd2V), as.double(ter_log),
       as.integer(model), as.double(A), as.double(W), as.integer(nq),
       as.double(eps), as.double(-1), as.double(rep(0, Ns)),
       PACKAGE = "diffIRT_orig"), times=30L)$time) / 1e6
  t_p <- median(microbenchmark(
    .C("LLdiff", rt_s, x_s, Ni, as.integer(nit),
       as.double(ai_log), as.double(vi),
       as.double(sd2A), as.double(sd2V), as.double(ter_log),
       as.integer(model), as.double(A), as.double(W), as.integer(nq),
       as.double(eps), as.double(-1), as.double(rep(0, Ns)),
       PACKAGE = "diffIRT_opt"), times=30L)$time) / 1e6
  cat(sprintf("N=%4d  orig=%.1fms  opt=%.1fms  speedup=%.2fx\n",
              Ns, t_o, t_p, t_o/t_p))
}
