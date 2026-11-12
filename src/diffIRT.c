/*
 *  diffIRT.c  -- optimized version
 *
 *  Changes from original (Dylan Molenaar, 2013):
 *
 *  1. LLdiff: the Navarro-Fuss series depends only on the boundary-separation
 *     node, so it is evaluated nq (not nq^2) times per response; the drift part
 *     is computed in log space. Log-sum-exp over quadrature nodes (no underflow).
 *  2. DERdiff: same forward-difference gradient as before, but per-item log
 *     densities are cached so each partial derivative recomputes only the
 *     item(s) the parameter touches, instead of the whole likelihood.
 *  3. OpenMP over subjects (deterministic serial reduction).
 *     Note: on macOS the default Apple toolchain has no OpenMP; the code then
 *     runs serially. src/Makevars uses $(SHLIB_OPENMP_CFLAGS).
 *  4. Hoisted loop invariants in LLfacA / LLfacV.
 *  5. factest: new FACscore (per-person 1-D Brent search, OpenMP) replaces the
 *     joint N-dimensional BFGS. The prior in LLfacA/LLfacV/FACscore is now a
 *     proper N(0, sd2) density; the original used exp(-.5*P/sqrt(sd2)).
 *  6. DERdiffAn: analytic gradient (series derivative computed with the
 *     same terms), selected from R by control$grad (default "analytic").
 *  7. Adaptive Gauss-Hermite quadrature (LLdiffAQ, DERdiffAQ, AQmoments):
 *     per-subject grids centred on the posterior; control$adaptive (default TRUE).
 *
 *  Original copyright 2013 Dylan Molenaar, University of Amsterdam.
 *  GNU GPL v2 or later.
 */

#include <math.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

/* ------------------------------------------------------------------ */
/* Forward declarations                                                 */
/* ------------------------------------------------------------------ */
void LLdiff(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL);
void DERdiff(double *y, double *x, int *N, int *nit, int *constrain,
             int *npar, double *ai, double *vi, double *sd2A, double *sd2V,
             double *ter, int *model, double *A, double *W, int *nq,
             double *epsAlt, double *delta, double *ou, double *LL,
             double *deriv);
void LLfacA(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, double *Ap, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL);
void LLfacV(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, double *Vp, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL);
void mdiffQQ(double *y, int *N, double *ai, double *vi, double *sd2A,
             double *sd2V, double *ter, int *model, double *A, double *W,
             int *nq, double *epsAlt, double *D);
void makeT(int *nit, double *design, int *dim, int *r, double *T);


/* ------------------------------------------------------------------ */
/* ddiff  (unchanged)                                                   */
/* ------------------------------------------------------------------ */
void ddiff(double *y, double *x, double *ap, double *ai,
           double *vp, double *vi, double *ter, int *model,
           double *epsAlt, double *ou)
{
  double aal = exp(*ap - *ai);
  double yy  = (*y - *ter) / (aal * aal);
  /* yy<=0 or NaN (y<=Ter, or a overflowed): density is 0. Without this guard,
     kl=Inf and (int)ceil(Inf) saturates to INT_MAX on ARM64 -> ~2e9 iterations. */
  if (!(yy > 0)) { *ou = 0.0; return; }
  double kl, ks, exl = 0.0;

  if (*model == 1) exl = (*vp - *vi);
  if (*model == 2) exl = exp(*vp - *vi);

  if ((M_PI * yy * (*epsAlt)) < 1) {
    kl = sqrt(-2.0 * log(M_PI * yy * (*epsAlt)) / (M_PI * M_PI * yy));
    if (kl < 1.0 / (M_PI * sqrt(yy))) kl = 1.0 / (M_PI * sqrt(yy));
  } else kl = 1.0 / (M_PI * sqrt(yy));

  if ((2.0 * sqrt(2.0 * M_PI * yy) * (*epsAlt)) < 1) {
    ks = 2.0 + sqrt(-2.0 * yy * log(2.0 * (*epsAlt) * sqrt(2.0 * M_PI * yy)));
    if (ks < (sqrt(yy) + 1.0)) ks = sqrt(yy) + 1.0;
  } else ks = 2.0;

  double p = 0.0;
  if (ks < kl) {
    for (int k = (int)(-floor((ks-1)/2)); k <= (int)ceil((ks-1)/2); k++)
      p += (.5 + 2*k) * exp(-(.5 + 2*k)*(.5 + 2*k) / (2.0*yy));
    p /= sqrt(2.0 * M_PI * yy * yy * yy);
  } else {
    int KK = (int)ceil(kl);
    for (int k = 1; k <= KK; k++)
      p += k * exp(-k*k * M_PI*M_PI * yy / 2.0) * sin(k * M_PI * .5);
    p *= M_PI;
  }
  p = p / (aal*aal) * exp(exl * aal * (*x - .5) - .5 * exl*exl * (*y - *ter));
  *ou = p;
}


/* ------------------------------------------------------------------ */
/* ddiffFacV  (unchanged)                                               */
/* ------------------------------------------------------------------ */
void ddiffFacV(double *y, double *x, double *ap, double *ai,
               double *vp, double *vi, double *ter, int *model,
               double *epsAlt, double *ou)
{
  double aal = exp(*ap - *ai);
  double yy  = (*y - *ter) / (aal * aal);
  /* yy<=0 or NaN (y<=Ter, or a overflowed): density is 0. Without this guard,
     kl=Inf and (int)ceil(Inf) saturates to INT_MAX on ARM64 -> ~2e9 iterations. */
  if (!(yy > 0)) { *ou = 0.0; return; }
  double kl, ks, exl = 0.0;

  if (*model == 1) exl = (*vp - *vi);
  if (*model == 2) exl = (*vp / exp(*vi));

  if ((M_PI * yy * (*epsAlt)) < 1) {
    kl = sqrt(-2.0 * log(M_PI * yy * (*epsAlt)) / (M_PI * M_PI * yy));
    if (kl < 1.0 / (M_PI * sqrt(yy))) kl = 1.0 / (M_PI * sqrt(yy));
  } else kl = 1.0 / (M_PI * sqrt(yy));

  if ((2.0 * sqrt(2.0 * M_PI * yy) * (*epsAlt)) < 1) {
    ks = 2.0 + sqrt(-2.0 * yy * log(2.0 * (*epsAlt) * sqrt(2.0 * M_PI * yy)));
    if (ks < (sqrt(yy) + 1.0)) ks = sqrt(yy) + 1.0;
  } else ks = 2.0;

  double p = 0.0;
  if (ks < kl) {
    for (int k = (int)(-floor((ks-1)/2)); k <= (int)ceil((ks-1)/2); k++)
      p += (.5 + 2*k) * exp(-(.5 + 2*k)*(.5 + 2*k) / (2.0*yy));
    p /= sqrt(2.0 * M_PI * yy * yy * yy);
  } else {
    int KK = (int)ceil(kl);
    for (int k = 1; k <= KK; k++)
      p += k * exp(-k*k * M_PI*M_PI * yy / 2.0) * sin(k * M_PI * .5);
    p *= M_PI;
  }
  p = p / (aal*aal) * exp(exl * aal * (*x - .5) - .5 * exl*exl * (*y - *ter));
  *ou = p;
}


/* ------------------------------------------------------------------ */
/* Fast marginal likelihood + gradient (LLdiff, DERdiff)                */
/*                                                                      */
/* Key facts exploited:                                                 */
/*  (a) The Navarro-Fuss series depends only on u = (y-Ter)/a^2, and    */
/*      a = exp(ap - ai) depends only on the boundary node n1, not on   */
/*      the drift node n2. So the series is evaluated nq (not nq^2)     */
/*      times per response; the drift part is cheap log-space algebra.  */
/*  (b) Everything is accumulated in log space, with log-sum-exp over   */
/*      quadrature nodes: no underflow for long tests.                  */
/*  (c) DERdiff keeps the original forward-difference gradient, but a   */
/*      perturbation of ai[j], vi[j] or ter[j] only changes item j's    */
/*      factor. Per-item log densities are cached, so each partial      */
/*      derivative costs one item (not a full likelihood). Perturbing   */
/*      vi[j] or sd2V does not even require re-evaluating the series.   */
/*  (d) OpenMP over subjects; per-subject results are summed serially,  */
/*      so results are deterministic regardless of thread count.        */
/* ------------------------------------------------------------------ */

#define LOGZERO (-1e10)          /* stand-in for log(0); keeps arithmetic finite */
#define LOGFLOOR (-690.7755279)  /* log(1e-300): the original per-subject floor  */

/* Navarro & Fuss (2009) series for the standardized density at w=.5,   */
/* u = t/a^2. Logic identical to the original ddiff().                  */
static double nf_series(double yy, double eps)
{
  double kl, ks;
  if ((M_PI * yy * eps) < 1) {
    kl = sqrt(-2.0 * log(M_PI * yy * eps) / (M_PI * M_PI * yy));
    if (kl < 1.0 / (M_PI * sqrt(yy))) kl = 1.0 / (M_PI * sqrt(yy));
  } else kl = 1.0 / (M_PI * sqrt(yy));
  if ((2.0 * sqrt(2.0 * M_PI * yy) * eps) < 1) {
    ks = 2.0 + sqrt(-2.0 * yy * log(2.0 * eps * sqrt(2.0 * M_PI * yy)));
    if (ks < (sqrt(yy) + 1.0)) ks = sqrt(yy) + 1.0;
  } else ks = 2.0;

  double p = 0.0;
  if (ks < kl) {
    for (int k = (int)(-floor((ks-1)/2)); k <= (int)ceil((ks-1)/2); k++)
      p += (.5 + 2*k) * exp(-(.5 + 2*k)*(.5 + 2*k) / (2.0*yy));
    p /= sqrt(2.0 * M_PI * yy * yy * yy);
  } else {
    int KK = (int)ceil(kl);
    for (int k = 1; k <= KK; k++)
      p += k * exp(-k*k * M_PI*M_PI * yy / 2.0) * sin(k * M_PI * .5);
    p *= M_PI;
  }
  return p;
}

/* log density of one response at all nq*nq nodes -> out[n1*nq + n2].  */
/* nodeA/nodeV are the scaled quadrature nodes for log(ap) and vp.      */
static void item_logf(double y, double x, double ai, double vi, double xter,
                      const double *nodeA, const double *nodeV, int nq,
                      int model, double eps, double *out)
{
  if (y == -999) {                       /* missing: contributes factor 1 */
    for (int q = 0; q < nq*nq; q++) out[q] = 0.0;
    return;
  }
  double t = y - xter, xm = x - .5;
  for (int n1 = 0; n1 < nq; n1++) {
    double la = nodeA[n1] - ai;          /* log a */
    double a  = exp(la);
    double yy = t / (a*a);
    double ls = LOGZERO;
    if (yy > 0) {
      double p = nf_series(yy, eps);
      if (p > 0 && R_FINITE(p)) ls = log(p) - 2.0*la;
    }
    double *o = out + n1*nq;
    if (ls == LOGZERO) { for (int n2 = 0; n2 < nq; n2++) o[n2] = LOGZERO; continue; }
    for (int n2 = 0; n2 < nq; n2++) {
      double v = (model == 1) ? (nodeV[n2] - vi) : exp(nodeV[n2] - vi);
      o[n2] = ls + v*a*xm - .5*v*v*t;
    }
  }
}

/* log marginal likelihood of one subject from per-node log products.   */
static double node_lse(const double *lp, const double *logw, int nq2)
{
  double m = -INFINITY;
  for (int q = 0; q < nq2; q++) { double z = lp[q] + logw[q]; if (z > m) m = z; }
  if (!(m > 0.5*LOGZERO) || !R_FINITE(m)) return LOGFLOOR;   /* all zero / NaN */
  double s = 0.0;
  for (int q = 0; q < nq2; q++) s += exp(lp[q] + logw[q] - m);
  double r = m + log(s);
  return R_FINITE(r) ? r : LOGFLOOR;
}

/* Per-subject: fill logf[j*nq2 + q] for all items and lp[q]; return logL. */
static double subject_eval(int i, const double *y, const double *x, int N, int nit,
                           const double *ai, const double *vi, const double *xter,
                           const double *nodeA, const double *nodeV, const double *logw,
                           int nq, int model, double eps, double *logf, double *lp)
{
  int nq2 = nq*nq;
  for (int q = 0; q < nq2; q++) lp[q] = 0.0;
  for (int j = 0; j < nit; j++) {
    double *f = logf + (size_t)j*nq2;
    item_logf(y[i + j*N], x[i + j*N], ai[j], vi[j], xter[j],
              nodeA, nodeV, nq, model, eps, f);
    for (int q = 0; q < nq2; q++) lp[q] += f[q];
  }
  return node_lse(lp, logw, nq2);
}

void LLdiff(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL)
{
  int n = *N, J = *nit, Q = *nq, nq2 = Q*Q;
  double *xter  = (double *)R_alloc(J, sizeof(double));
  double *nodeA = (double *)R_alloc(Q, sizeof(double));
  double *nodeV = (double *)R_alloc(Q, sizeof(double));
  double *logw  = (double *)R_alloc(nq2, sizeof(double));
  double *logL  = (double *)R_alloc(n, sizeof(double));
  for (int j = 0; j < J; j++) xter[j] = exp(ter[j]);
  for (int q = 0; q < Q; q++) { nodeA[q] = sqrt(exp(*sd2A))*A[q]; nodeV[q] = sqrt(exp(*sd2V))*A[q]; }
  for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) logw[a*Q+b] = log(W[a]) + log(W[b]);

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    double *logf = (double *)malloc((size_t)J*nq2*sizeof(double));
    double *lp   = (double *)malloc((size_t)nq2*sizeof(double));
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int i = 0; i < n; i++)
      logL[i] = subject_eval(i, y, x, n, J, ai, vi, xter, nodeA, nodeV, logw,
                             Q, *model, *epsAlt, logf, lp);
    free(logf); free(lp);
  }

  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += logL[i];
    LL[i] = (logL[i] <= LOGFLOOR) ? 1e-300 : exp(logL[i]);
  }
  *ou = -2.0 * s;
}

void DERdiff(double *y, double *x, int *N, int *nit, int *constrain,
             int *npar, double *ai, double *vi, double *sd2A, double *sd2V,
             double *ter, int *model, double *A, double *W, int *nq,
             double *epsAlt, double *delta, double *ou, double *LL,
             double *deriv)
{
  int n = *N, J = *nit, Q = *nq, nq2 = Q*Q, P = *npar;
  double d = *delta;
  double *xter   = (double *)R_alloc(J, sizeof(double));
  double *xterP  = (double *)R_alloc(J, sizeof(double));   /* exp(ter + delta) */
  double *nodeA  = (double *)R_alloc(Q, sizeof(double));
  double *nodeV  = (double *)R_alloc(Q, sizeof(double));
  double *nodeAP = (double *)R_alloc(Q, sizeof(double));   /* sd2A + delta */
  double *nodeVP = (double *)R_alloc(Q, sizeof(double));   /* sd2V + delta */
  double *logw   = (double *)R_alloc(nq2, sizeof(double));
  double *logL   = (double *)R_alloc(n, sizeof(double));
  double *dL     = (double *)R_alloc((size_t)n*P, sizeof(double));
  for (int j = 0; j < J; j++) { xter[j] = exp(ter[j]); xterP[j] = exp(ter[j] + d); }
  for (int q = 0; q < Q; q++) {
    nodeA[q]  = sqrt(exp(*sd2A))*A[q];     nodeV[q]  = sqrt(exp(*sd2V))*A[q];
    nodeAP[q] = sqrt(exp(*sd2A + d))*A[q]; nodeVP[q] = sqrt(exp(*sd2V + d))*A[q];
  }
  for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) logw[a*Q+b] = log(W[a]) + log(W[b]);
  const int *cA = constrain, *cV = constrain + J, *cT = constrain + 2*J;
  int csA = constrain[3*J], csV = constrain[3*J+1];

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    double *logf = (double *)malloc((size_t)J*nq2*sizeof(double));
    double *lp   = (double *)malloc((size_t)nq2*sizeof(double));
    double *lpP  = (double *)malloc((size_t)nq2*sizeof(double));
    double *fP   = (double *)malloc((size_t)nq2*sizeof(double));
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int i = 0; i < n; i++) {
      double l0 = subject_eval(i, y, x, n, J, ai, vi, xter, nodeA, nodeV, logw,
                               Q, *model, *epsAlt, logf, lp);
      logL[i] = l0;
      for (int p = 1; p <= P; p++) {
        int pA = (csA == p), pV = (csV == p), any = 0;
        const double *nA = pA ? nodeAP : nodeA, *nV = pV ? nodeVP : nodeV;
        for (int q = 0; q < nq2; q++) lpP[q] = lp[q];
        for (int j = 0; j < J; j++) {
          int ja = (cA[j] == p), jv = (cV[j] == p), jt = (cT[j] == p);
          if (!(ja || jv || jt || pA || pV)) continue;
          any = 1;
          const double *f = logf + (size_t)j*nq2;
          item_logf(y[i + j*n], x[i + j*n], ai[j] + (ja ? d : 0.0), vi[j] + (jv ? d : 0.0),
                    jt ? xterP[j] : xter[j], nA, nV, Q, *model, *epsAlt, fP);
          for (int q = 0; q < nq2; q++) lpP[q] += fP[q] - f[q];
        }
        dL[(size_t)i*P + (p-1)] = any ? -2.0*(node_lse(lpP, logw, nq2) - l0)/d : 0.0;
      }
    }
    free(logf); free(lp); free(lpP); free(fP);
  }

  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += logL[i];
    LL[i] = (logL[i] <= LOGFLOOR) ? 1e-300 : exp(logL[i]);
  }
  *ou = -2.0 * s;
  for (int p = 0; p < P; p++) {
    double g = 0.0;
    for (int i = 0; i < n; i++) g += dL[(size_t)i*P + p];
    deriv[p] = g;
  }
}

/* ------------------------------------------------------------------ */
/* Analytic gradient                                                    */
/* ------------------------------------------------------------------ */

/* nf_series plus d log p / du, using the same terms (so the derivative   */
/* is exact for the truncated series between its term-count jumps).      */
static double nf_series_d(double yy, double eps, double *dlogp)
{
  double kl, ks;
  if ((M_PI * yy * eps) < 1) {
    kl = sqrt(-2.0 * log(M_PI * yy * eps) / (M_PI * M_PI * yy));
    if (kl < 1.0 / (M_PI * sqrt(yy))) kl = 1.0 / (M_PI * sqrt(yy));
  } else kl = 1.0 / (M_PI * sqrt(yy));
  if ((2.0 * sqrt(2.0 * M_PI * yy) * eps) < 1) {
    ks = 2.0 + sqrt(-2.0 * yy * log(2.0 * eps * sqrt(2.0 * M_PI * yy)));
    if (ks < (sqrt(yy) + 1.0)) ks = sqrt(yy) + 1.0;
  } else ks = 2.0;

  double p = 0.0, dp = 0.0;
  if (ks < kl) {
    /* p = S / sqrt(2 pi u^3),  S = sum c exp(-c^2 / 2u) */
    double S = 0.0, dS = 0.0;
    for (int k = (int)(-floor((ks-1)/2)); k <= (int)ceil((ks-1)/2); k++) {
      double c = .5 + 2*k, e = c * exp(-c*c / (2.0*yy));
      S += e; dS += e * c*c / (2.0*yy*yy);
    }
    double nrm = sqrt(2.0 * M_PI * yy * yy * yy);
    p = S / nrm;
    dp = dS / nrm - 1.5 * p / yy;
  } else {
    int KK = (int)ceil(kl);
    for (int k = 1; k <= KK; k++) {
      double e = k * exp(-k*k * M_PI*M_PI * yy / 2.0) * sin(k * M_PI * .5);
      p += e; dp += e * (-k*k * M_PI*M_PI / 2.0);
    }
    p *= M_PI; dp *= M_PI;
  }
  *dlogp = dp / p;
  return p;
}

/* Gradient of -2 logL with respect to the full parameter vector          */
/* (ai[1..J], vi[1..J], ter[1..J], sd2A, sd2V), summed into `deriv` by     */
/* `constrain` (fixed parameters, constrain == 0, are skipped). Same       */
/* signature as DERdiff; `delta` is unused.                                */
/*                                                                       */
/* With posterior node weights w_q = W_q f_q / sum(W f), d logL_i / dθ =   */
/* sum_q w_q sum_j d log f_ij(q) / dθ. For one response,                   */
/*   log f = log p(u) - 2 la + v a xm - .5 v^2 t,  u = t / a^2, la = log a  */
/*   d/dla = -2 u p'/p - 2 + v a xm,  d/dt = p'/p / a^2 - .5 v^2,          */
/*   d/dv  = a xm - v t.                                                   */
void DERdiffAn(double *y, double *x, int *N, int *nit, int *constrain,
               int *npar, double *ai, double *vi, double *sd2A, double *sd2V,
               double *ter, int *model, double *A, double *W, int *nq,
               double *epsAlt, double *delta, double *ou, double *LL,
               double *deriv)
{
  int n = *N, J = *nit, Q = *nq, nq2 = Q*Q, P = *npar, K = 3*J + 2, mod = *model;
  double eps = *epsAlt;
  double *xter  = (double *)R_alloc(J, sizeof(double));
  double *nodeA = (double *)R_alloc(Q, sizeof(double));
  double *nodeV = (double *)R_alloc(Q, sizeof(double));
  double *logw  = (double *)R_alloc(nq2, sizeof(double));
  double *logL  = (double *)R_alloc(n, sizeof(double));
  double *G     = (double *)R_alloc((size_t)n*K, sizeof(double));   /* per-subject d logL */
  for (int j = 0; j < J; j++) xter[j] = exp(ter[j]);
  for (int q = 0; q < Q; q++) { nodeA[q] = sqrt(exp(*sd2A))*A[q]; nodeV[q] = sqrt(exp(*sd2V))*A[q]; }
  for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) logw[a*Q+b] = log(W[a]) + log(W[b]);

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    /* per (item, boundary node): log p - 2 la, d log p/du, a, u */
    double *ls = (double *)malloc((size_t)J*Q*sizeof(double));
    double *dl = (double *)malloc((size_t)J*Q*sizeof(double));
    double *aa = (double *)malloc((size_t)J*Q*sizeof(double));
    double *uu = (double *)malloc((size_t)J*Q*sizeof(double));
    double *lp = (double *)malloc((size_t)nq2*sizeof(double));
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int i = 0; i < n; i++) {
      double *g = G + (size_t)i*K;
      for (int k = 0; k < K; k++) g[k] = 0.0;
      for (int q = 0; q < nq2; q++) lp[q] = 0.0;

      /* pass 1: log density at every node (as item_logf), caching pieces */
      for (int j = 0; j < J; j++) {
        double yj = y[i + j*n];
        if (yj == -999) continue;
        double t = yj - xter[j], xm = x[i + j*n] - .5;
        for (int n1 = 0; n1 < Q; n1++) {
          int c = j*Q + n1;
          double la = nodeA[n1] - ai[j], a = exp(la), u = t/(a*a), d = 0.0;
          ls[c] = LOGZERO; aa[c] = a; uu[c] = u;
          if (u > 0) {
            double p = nf_series_d(u, eps, &d);
            if (p > 0 && R_FINITE(p) && R_FINITE(d)) ls[c] = log(p) - 2.0*la;
          }
          dl[c] = d;
          double *l = lp + n1*Q;
          if (ls[c] == LOGZERO) { for (int n2 = 0; n2 < Q; n2++) l[n2] += LOGZERO; continue; }
          for (int n2 = 0; n2 < Q; n2++) {
            double v = (mod == 1) ? (nodeV[n2] - vi[j]) : exp(nodeV[n2] - vi[j]);
            l[n2] += ls[c] + v*a*xm - .5*v*v*t;
          }
        }
      }
      double l0 = node_lse(lp, logw, nq2);
      logL[i] = l0;
      if (l0 <= LOGFLOOR) continue;          /* floored subject: flat, gradient 0 */
      for (int q = 0; q < nq2; q++) lp[q] = exp(lp[q] + logw[q] - l0);   /* posterior weights */

      /* pass 2: accumulate weighted derivatives */
      double gsA = 0.0, gsV = 0.0;
      for (int j = 0; j < J; j++) {
        double yj = y[i + j*n];
        if (yj == -999) continue;
        double t = yj - xter[j], xm = x[i + j*n] - .5;
        double gai = 0.0, gvi = 0.0, gt = 0.0;
        for (int n1 = 0; n1 < Q; n1++) {
          int c = j*Q + n1;
          if (ls[c] == LOGZERO) continue;    /* node weight is 0 */
          double a = aa[c], u = uu[c], d = dl[c];
          double dla0 = -2.0*u*d - 2.0, dt0 = d/(a*a);
          const double *w = lp + n1*Q;
          double sla = 0.0;
          for (int n2 = 0; n2 < Q; n2++) {
            double wq = w[n2];
            if (wq == 0.0) continue;
            double v = (mod == 1) ? (nodeV[n2] - vi[j]) : exp(nodeV[n2] - vi[j]);
            double dla = dla0 + v*a*xm, dv = a*xm - v*t, dt = dt0 - .5*v*v;
            double dvdvi = (mod == 1) ? -1.0 : -v;           /* v = node - vi  or  exp(node - vi) */
            double dvdsV = (mod == 1) ? .5*nodeV[n2] : .5*v*nodeV[n2];
            sla += wq*dla;
            gvi += wq*dv*dvdvi;
            gt  += wq*dt;
            gsV += wq*dv*dvdsV;
          }
          gai -= sla;                         /* la = nodeA - ai */
          gsA += sla*.5*nodeA[n1];            /* nodeA = exp(sd2A/2) A */
        }
        g[j] = gai; g[J + j] = gvi; g[2*J + j] = -gt*xter[j];   /* t = y - exp(ter) */
      }
      g[3*J] = gsA; g[3*J + 1] = gsV;
    }
    free(ls); free(dl); free(aa); free(uu); free(lp);
  }

  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += logL[i];
    LL[i] = (logL[i] <= LOGFLOOR) ? 1e-300 : exp(logL[i]);
  }
  *ou = -2.0 * s;
  for (int p = 0; p < P; p++) deriv[p] = 0.0;
  for (int k = 0; k < K; k++) {
    int p = constrain[k];
    if (p < 1 || p > P) continue;
    double gk = 0.0;
    for (int i = 0; i < n; i++) gk += G[(size_t)i*K + k];
    deriv[p-1] += -2.0*gk;
  }
}

/* ------------------------------------------------------------------ */
/* Adaptive Gauss-Hermite quadrature                                    */
/* ------------------------------------------------------------------ */
/* Subject i gets its own product grid: log-boundary nodes               */
/*   zA = muA[i] + sA[i] * A[n1],  drift nodes zV = muV[i] + sV[i] * A[n2] */
/* (A, W: Gauss-Hermite rule for a standard normal). With prior           */
/* N(0, s2A) x N(0, s2V), the marginal likelihood is                       */
/*   sum_q W[n1] W[n2] * prior(zA, zV) / proposal(zA, zV) * prod_j f_ij    */
/*   log(prior/proposal) for one axis = log(s/sigma) - z^2/(2 sigma^2) + A^2/2. */
/* With mu = 0 and s = sigma this is exactly the non-adaptive rule. The    */
/* nodes are held fixed while the optimizer runs (so the objective is a    */
/* fixed smooth function and the gradient is exact); the population        */
/* variances then enter only through the prior term.                      */
/* A posterior SD never exceeds the prior SD, so each grid scale is capped */
/* at the current prior SD. Without the cap the rule is invalid once       */
/* sigma < s (the weight s/sigma at the centre node grows without bound as */
/* sigma -> 0, which the optimizer can exploit). Where the cap is active   */
/* the nodes move with sigma, and the gradient includes that dependence.   */

/* Workspace per thread: 4*J*Q + Q*Q + 6*Q doubles. Returns log L_i.      */
/* g (length 3J+2, or NULL): d log L_i / d (ai, vi, ter, sd2A, sd2V).      */
/* mom (length 4, or NULL): posterior mean and SD of log ap and of vp.    */
static double aq_subject(int i, const double *y, const double *x, int n, int J,
                         const double *ai, const double *vi, const double *xter,
                         double s2A, double s2V, int mod,
                         const double *A, const double *W, int Q, double eps,
                         double muA, double sA, double muV, double sV,
                         double *ws, double *g, double *mom)
{
  int nq2 = Q*Q;
  double *ls = ws, *dl = ls + J*Q, *aa = dl + J*Q, *uu = aa + J*Q, *lp = uu + J*Q;
  double *zA = lp + nq2, *zV = zA + Q, *lwA = zV + Q, *lwV = lwA + Q;
  double *accA = lwV + Q, *accV = accA + Q;
  double sigA = sqrt(s2A), sigV = sqrt(s2V);
  /* a variance fixed at (or underflowing to) 0: all nodes at 0, as in the */
  /* non-adaptive rule, and no prior term                                  */
  int degA = !(s2A > 0), degV = !(s2V > 0);
  int capA = !degA && (sA > sigA), capV = !degV && (sV > sigV);
  if (capA) sA = sigA;
  if (capV) sV = sigV;
  for (int k = 0; k < Q; k++) {
    zA[k] = degA ? 0.0 : muA + sA*A[k];
    zV[k] = degV ? 0.0 : muV + sV*A[k];
    lwA[k] = degA ? log(W[k]) : log(W[k]) + log(sA) - .5*log(s2A) - zA[k]*zA[k]/(2.0*s2A) + .5*A[k]*A[k];
    lwV[k] = degV ? log(W[k]) : log(W[k]) + log(sV) - .5*log(s2V) - zV[k]*zV[k]/(2.0*s2V) + .5*A[k]*A[k];
  }
  for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) lp[a*Q + b] = lwA[a] + lwV[b];

  /* pass 1: log density of every response at every node */
  for (int j = 0; j < J; j++) {
    double yj = y[i + j*n];
    if (yj == -999) continue;
    double t = yj - xter[j], xm = x[i + j*n] - .5;
    for (int n1 = 0; n1 < Q; n1++) {
      int c = j*Q + n1;
      double la = zA[n1] - ai[j], a = exp(la), u = t/(a*a), d = 0.0;
      ls[c] = LOGZERO; aa[c] = a; uu[c] = u;
      if (u > 0) {
        double p = g ? nf_series_d(u, eps, &d) : nf_series(u, eps);
        if (p > 0 && R_FINITE(p) && R_FINITE(d)) ls[c] = log(p) - 2.0*la;
      }
      dl[c] = d;
      double *l = lp + n1*Q;
      if (ls[c] == LOGZERO) { for (int n2 = 0; n2 < Q; n2++) l[n2] += LOGZERO; continue; }
      for (int n2 = 0; n2 < Q; n2++) {
        double v = (mod == 1) ? (zV[n2] - vi[j]) : exp(zV[n2] - vi[j]);
        l[n2] += ls[c] + v*a*xm - .5*v*v*t;
      }
    }
  }
  double m = -INFINITY;
  for (int q = 0; q < nq2; q++) if (lp[q] > m) m = lp[q];
  double l0 = LOGFLOOR;
  if (m > 0.5*LOGZERO && R_FINITE(m)) {
    double s = 0.0;
    for (int q = 0; q < nq2; q++) s += exp(lp[q] - m);
    l0 = m + log(s);
    if (!R_FINITE(l0)) l0 = LOGFLOOR;
  }
  if (g) for (int k = 0; k < 3*J + 2; k++) g[k] = 0.0;
  if (l0 <= LOGFLOOR) {                    /* floored subject: flat, keep old nodes */
    if (mom) { mom[0] = muA; mom[1] = sA; mom[2] = muV; mom[3] = sV; }
    return LOGFLOOR;
  }
  for (int q = 0; q < nq2; q++) lp[q] = exp(lp[q] - l0);   /* posterior weights */

  if (mom) {
    double mA = 0, mV = 0, vA = 0, vV = 0;
    for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) {
      double w = lp[a*Q + b]; mA += w*zA[a]; mV += w*zV[b];
    }
    for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) {
      double w = lp[a*Q + b]; vA += w*(zA[a]-mA)*(zA[a]-mA); vV += w*(zV[b]-mV)*(zV[b]-mV);
    }
    mom[0] = mA; mom[1] = sqrt(vA); mom[2] = mV; mom[3] = sqrt(vV);
  }
  if (!g) return l0;

  /* pass 2: gradient (items as in DERdiffAn; variances via the prior term) */
  /* d log(prior/proposal weight) / d log sigma^2. Uncapped: nodes fixed,    */
  /* -1/2 + z^2/(2 s2). Capped (s = sigma, z = mu + sigma A): the log(s/sigma) */
  /* term is constant and dz/d log s2 = sigma A / 2, giving                 */
  /* -z A/(2 sigma) + z^2/(2 s2); the items' dependence on z is added below. */
  double gsA = 0.0, gsV = 0.0;
  for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) {
    double w = lp[a*Q + b];
    if (!degA) gsA += w*((capA ? -zA[a]*A[a]/(2.0*sigA) : -.5) + zA[a]*zA[a]/(2.0*s2A));
    if (!degV) gsV += w*((capV ? -zV[b]*A[b]/(2.0*sigV) : -.5) + zV[b]*zV[b]/(2.0*s2V));
  }
  for (int k = 0; k < Q; k++) { accA[k] = 0.0; accV[k] = 0.0; }
  for (int j = 0; j < J; j++) {
    double yj = y[i + j*n];
    if (yj == -999) continue;
    double t = yj - xter[j], xm = x[i + j*n] - .5;
    double gai = 0.0, gvi = 0.0, gt = 0.0;
    for (int n1 = 0; n1 < Q; n1++) {
      int c = j*Q + n1;
      if (ls[c] == LOGZERO) continue;
      double a = aa[c], u = uu[c], d = dl[c];
      double dla0 = -2.0*u*d - 2.0, dt0 = d/(a*a);
      const double *w = lp + n1*Q;
      for (int n2 = 0; n2 < Q; n2++) {
        double wq = w[n2];
        if (wq == 0.0) continue;
        double v = (mod == 1) ? (zV[n2] - vi[j]) : exp(zV[n2] - vi[j]);
        double dv = a*xm - v*t, dvdz = (mod == 1) ? 1.0 : v;   /* v = zV - vi or exp(zV - vi) */
        double dla = dla0 + v*a*xm;                            /* = d log f / d zA */
        gai -= wq*dla;
        gvi += wq*dv*(-dvdz);
        gt  += wq*(dt0 - .5*v*v);
        accA[n1] += wq*dla;
        accV[n2] += wq*dv*dvdz;
      }
    }
    g[j] = gai; g[J + j] = gvi; g[2*J + j] = -gt*xter[j];
  }
  if (capA) for (int k = 0; k < Q; k++) gsA += accA[k]*A[k]*.5*sigA;
  if (capV) for (int k = 0; k < Q; k++) gsV += accV[k]*A[k]*.5*sigV;
  g[3*J] = gsA; g[3*J + 1] = gsV;
  return l0;
}

/* Shared driver. what: 0 = likelihood, 1 = gradient, 2 = posterior moments. */
static void aq_driver(int what, double *y, double *x, int *N, int *nit, int *constrain, int *npar,
                      double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model,
                      double *A, double *W, int *nq, double *epsAlt,
                      double *muA, double *sA, double *muV, double *sV,
                      double *ou, double *LL, double *deriv, double *mom)
{
  int n = *N, J = *nit, Q = *nq, K = 3*J + 2;
  double s2A = exp(*sd2A), s2V = exp(*sd2V);
  double *xter = (double *)R_alloc(J, sizeof(double));
  double *logL = (double *)R_alloc(n, sizeof(double));
  double *G = (what == 1) ? (double *)R_alloc((size_t)n*K, sizeof(double)) : NULL;
  for (int j = 0; j < J; j++) xter[j] = exp(ter[j]);
  size_t wsz = (size_t)4*J*Q + (size_t)Q*Q + 6*(size_t)Q;

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    double *ws = (double *)malloc(wsz*sizeof(double));
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int i = 0; i < n; i++)
      logL[i] = aq_subject(i, y, x, n, J, ai, vi, xter, s2A, s2V, *model, A, W, Q, *epsAlt,
                           muA[i], sA[i], muV[i], sV[i], ws,
                           G ? G + (size_t)i*K : NULL, (what == 2) ? mom + 4*(size_t)i : NULL);
    free(ws);
  }

  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += logL[i];
    if (LL) LL[i] = (logL[i] <= LOGFLOOR) ? 1e-300 : exp(logL[i]);
  }
  if (ou) *ou = -2.0 * s;
  if (what == 1) {
    int P = *npar;
    for (int p = 0; p < P; p++) deriv[p] = 0.0;
    for (int k = 0; k < K; k++) {
      int p = constrain[k];
      if (p < 1 || p > P) continue;
      double gk = 0.0;
      for (int i = 0; i < n; i++) gk += G[(size_t)i*K + k];
      deriv[p-1] += -2.0*gk;
    }
  }
}

void LLdiffAQ(double *y, double *x, int *N, int *nit, double *ai, double *vi,
              double *sd2A, double *sd2V, double *ter, int *model,
              double *A, double *W, int *nq, double *epsAlt,
              double *muA, double *sA, double *muV, double *sV, double *ou, double *LL)
{
  aq_driver(0, y, x, N, nit, NULL, NULL, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt,
            muA, sA, muV, sV, ou, LL, NULL, NULL);
}

void DERdiffAQ(double *y, double *x, int *N, int *nit, int *constrain, int *npar,
               double *ai, double *vi, double *sd2A, double *sd2V, double *ter, int *model,
               double *A, double *W, int *nq, double *epsAlt,
               double *muA, double *sA, double *muV, double *sV,
               double *ou, double *LL, double *deriv)
{
  aq_driver(1, y, x, N, nit, constrain, npar, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt,
            muA, sA, muV, sV, ou, LL, deriv, NULL);
}

/* Posterior moments (mom[4*i + 0..3] = mean, SD of log ap; mean, SD of vp) */
/* computed on the current grid; used to re-centre the grid.               */
void AQmoments(double *y, double *x, int *N, int *nit, double *ai, double *vi,
               double *sd2A, double *sd2V, double *ter, int *model,
               double *A, double *W, int *nq, double *epsAlt,
               double *muA, double *sA, double *muV, double *sV, double *ou, double *mom)
{
  aq_driver(2, y, x, N, nit, NULL, NULL, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt,
            muA, sA, muV, sV, ou, NULL, NULL, mom);
}

/* ------------------------------------------------------------------ */
/* LLfacA  -- OPTIMIZED (hoisted + log-space + OpenMP)                 */
/* ------------------------------------------------------------------ */
void LLfacA(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, double *Ap, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL)
{
  double sqrtSdV = sqrt(exp(*sd2V));

  double *xter_cache = (double *)R_alloc(*nit, sizeof(double));
  for (int j = 0; j < *nit; j++) xter_cache[j] = exp(ter[j]);

  double *AtransV = (double *)R_alloc(*nq, sizeof(double));
  for (int n = 0; n < *nq; n++) AtransV[n] = sqrtSdV * A[n];

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 4)
  #endif
  for (int i = 0; i < *N; i++) {
    double tmp2 = 0.0;

    for (int n1 = 0; n1 < *nq; n1++) {
      double aVp = AtransV[n1];
      double log_tmp1 = 0.0;
      int valid = 1;

      for (int j = 0; j < *nit; j++) {
        if (y[i + j*(*N)] != -999) {
          double ou_local;
          ddiff(&y[i + j*(*N)], &x[i + j*(*N)],
                &Ap[i], &ai[j], &aVp, &vi[j],
                &xter_cache[j], model, epsAlt, &ou_local);
          if (ou_local <= 0.0) { valid = 0; break; }
          log_tmp1 += log(ou_local);
        }
      }
      if (valid) tmp2 += exp(log_tmp1) * W[n1];
    }

    if (tmp2 == 0.0 || !R_FINITE(tmp2))
      LL[i] = 1e-300;
    else
      LL[i] = tmp2 / sqrt(2.0 * M_PI * exp(*sd2A))
               * exp(-0.5 * Ap[i] * Ap[i] / exp(*sd2A));   /* N(0, sd2A) prior on log a_p */
  }

  double logsum = 0.0;
  for (int i = 0; i < *N; i++) logsum += log(LL[i]);
  *ou = -2.0 * logsum;
}


/* ------------------------------------------------------------------ */
/* LLfacV  -- OPTIMIZED (hoisted + log-space + OpenMP)                 */
/* ------------------------------------------------------------------ */
void LLfacV(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, double *Vp, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL)
{
  double sqrtSdA = sqrt(exp(*sd2A));

  double *xter_cache = (double *)R_alloc(*nit, sizeof(double));
  for (int j = 0; j < *nit; j++) xter_cache[j] = exp(ter[j]);

  double *AtransA = (double *)R_alloc(*nq, sizeof(double));
  for (int n = 0; n < *nq; n++) AtransA[n] = sqrtSdA * A[n];

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 4)
  #endif
  for (int i = 0; i < *N; i++) {
    double tmp2 = 0.0;

    for (int n1 = 0; n1 < *nq; n1++) {
      double aAp = AtransA[n1];
      double log_tmp1 = 0.0;
      int valid = 1;

      for (int j = 0; j < *nit; j++) {
        if (y[i + j*(*N)] != -999) {
          double ou_local;
          ddiffFacV(&y[i + j*(*N)], &x[i + j*(*N)],
                    &aAp, &ai[j], &Vp[i], &vi[j],
                    &xter_cache[j], model, epsAlt, &ou_local);
          if (ou_local <= 0.0) { valid = 0; break; }
          log_tmp1 += log(ou_local);
        }
      }
      if (valid) tmp2 += exp(log_tmp1) * W[n1];
    }

    if (tmp2 == 0.0 || !R_FINITE(tmp2))
      LL[i] = 1e-300;
    else {
      double PV = (*model == 1) ? Vp[i] : log(Vp[i]);
      LL[i] = tmp2 / sqrt(2.0 * M_PI * exp(*sd2V))
               * exp(-0.5 * PV * PV / exp(*sd2V));          /* N(0, sd2V) prior on v_p (D) or log theta_p (Q) */
    }
  }

  double logsum = 0.0;
  for (int i = 0; i < *N; i++) logsum += log(LL[i]);
  *ou = -2.0 * logsum;
}


/* ------------------------------------------------------------------ */
/* ddiffQQ  (unchanged)                                                 */
/* ------------------------------------------------------------------ */
void ddiffQQ(double *y, double *ap, double *ai, double *vp, double *vi,
             double *ter, int *model, double *epsAlt, double *ou)
{
  double aal = exp(*ap - *ai);
  double yy  = (*y - *ter) / (aal * aal);
  /* yy<=0 or NaN (y<=Ter, or a overflowed): density is 0. Without this guard,
     kl=Inf and (int)ceil(Inf) saturates to INT_MAX on ARM64 -> ~2e9 iterations. */
  if (!(yy > 0)) { *ou = 0.0; return; }
  double kl, ks, exl = 0.0;

  if (*model == 1) exl = (*vp - *vi);
  if (*model == 2) exl = exp(*vp - *vi);

  if ((M_PI * yy * (*epsAlt)) < 1) {
    kl = sqrt(-2.0 * log(M_PI * yy * (*epsAlt)) / (M_PI * M_PI * yy));
    if (kl < 1.0 / (M_PI * sqrt(yy))) kl = 1.0 / (M_PI * sqrt(yy));
  } else kl = 1.0 / (M_PI * sqrt(yy));

  if ((2.0 * sqrt(2.0 * M_PI * yy) * (*epsAlt)) < 1) {
    ks = 2.0 + sqrt(-2.0 * yy * log(2.0 * (*epsAlt) * sqrt(2.0 * M_PI * yy)));
    if (ks < (sqrt(yy) + 1.0)) ks = sqrt(yy) + 1.0;
  } else ks = 2.0;

  double p = 0.0;
  if (ks < kl) {
    for (int k = (int)(-floor((ks-1)/2)); k <= (int)ceil((ks-1)/2); k++)
      p += (.5 + 2*k) * exp(-(.5 + 2*k)*(.5 + 2*k) / (2.0*yy));
    p /= sqrt(2.0 * M_PI * yy * yy * yy);
  } else {
    int KK = (int)ceil(kl);
    for (int k = 1; k <= KK; k++)
      p += k * exp(-k*k * M_PI*M_PI * yy / 2.0) * sin(k * M_PI * .5);
    p *= M_PI;
  }
  /* density of both responses, in log space: exp(exl*aal/2) alone overflows */
  /* for large drift and boundary, although the product with exp(-exl^2 t/2) */
  /* is small                                                                */
  if (!(p > 0)) { *ou = 0.0; return; }
  double lb = log(p) - 2.0*log(aal) - .5*exl*exl*(*y-*ter);
  *ou = exp(lb - .5*exl*aal) + exp(lb + .5*exl*aal);
}


/* ------------------------------------------------------------------ */
/* mdiffQQ: marginal RT density of one item (both responses)           */
/* ------------------------------------------------------------------ */
void mdiffQQ(double *y, int *N, double *ai, double *vi, double *sd2A,
             double *sd2V, double *ter, int *model, double *A, double *W,
             int *nq, double *epsAlt, double *D)
{
  double xter, tmp, A1trans, A2trans;
  double ou = 0.0;

  for (int i = 0; i < *N; i++) {
    tmp = 0.0;
    for (int n1 = 0; n1 < *nq; n1++) for (int n2 = 0; n2 < *nq; n2++) {
      A1trans = sqrt(exp(*sd2A)) * A[n1];
      A2trans = sqrt(exp(*sd2V)) * A[n2];
      if (y[i] != -999) {
        xter = exp(*ter);
        ddiffQQ(&y[i], &A1trans, ai, &A2trans, vi, &xter, model, epsAlt, &ou);
      }
      tmp += ou * W[n1] * W[n2];
    }
    /* (was -999 for a zero density, which pdiffIRT then integrated as a value) */
    D[i] = R_FINITE(tmp) ? tmp : NA_REAL;
  }
}


/* ------------------------------------------------------------------ */
/* equal, makeT  (unchanged)                                            */
/* ------------------------------------------------------------------ */
void equal(double *big, double *small, int *dim, int *ou)
{
  int d = 0;
  int lS = (int)small[0];
  int lB = (int)big[0];
  for (int s = 1; s <= lS; s++)
    for (int b = 1; b <= lB; b++)
      if (big[*dim * b] == small[*dim * s]) d++;
  *ou = (d == lS) ? 1 : 0;
}

void makeT(int *nit, double *design, int *dim, int *r, double *T)
{
  int ou;
  for (int i = 0; i < *r; i++)
    for (int j = 0; j < *dim; j++)
      if (i <= j) {
        equal(&design[j], &design[i], dim, &ou);
        T[i + j*(*r)] = ou;
      }
}


/* ------------------------------------------------------------------ */
/* FACscore: per-subject factor scores for factest()                    */
/*                                                                      */
/* The old factest() ran one BFGS over all N person parameters with a   */
/* numerical gradient (N+1 full LLfacA/LLfacV calls per gradient, and   */
/* an N x N BFGS matrix). The objective is a sum of independent         */
/* per-subject terms, so here each subject is minimized on its own      */
/* with Brent's method (OpenMP over subjects).                          */
/*                                                                      */
/* Per-subject objective is identical to the corresponding term of      */
/* LLfacA (which = 1, gamma) or LLfacV (which = 2, theta), including    */
/* its prior term.                                                      */
/*   which = 1: par = log gamma_p (Ap); integrate over drift nodes.     */
/*   which = 2: par = Vp (model D) or log theta_p (model Q, Vp=exp(par))*/
/*              integrate over boundary nodes.                          */
/* Outputs: est[i], hess[i] (numerical 2nd derivative of -2 log post),  */
/* conv[i] (0 = interior optimum, 1 = at search-interval edge).         */
/* adapt = 1: integrate over the other latent variable with each        */
/* subject's adaptive grid (muO[i] + sO[i] * A, with prior/proposal     */
/* weights and the scale capped at the prior SD, as in aq_subject);     */
/* adapt = 0: the non-adaptive grid sd * A.                             */
/* ------------------------------------------------------------------ */

typedef struct {
  const double *y, *x, *ai, *vi, *xter, *node, *logw;
  int N, nit, nq, model, which, i;
  double sd2A, sd2V, eps;
  double *ls, *a, *lp;     /* workspace: ls/a are nit*nq (which = 2 cache) */
} facinfo;

static double fac_obj(double par, void *vinfo)
{
  facinfo *f = (facinfo *)vinfo;
  int J = f->nit, Q = f->nq, N = f->N, i = f->i;
  for (int n = 0; n < Q; n++) f->lp[n] = 0.0;
  double s2, Pval;
  if (f->which == 1) {
    double ap = par;
    for (int j = 0; j < J; j++) {
      double yj = f->y[i + j*N];
      if (yj == -999) continue;
      double t = yj - f->xter[j], xm = f->x[i + j*N] - .5;
      double la = ap - f->ai[j], a = exp(la), yy = t/(a*a), ls = LOGZERO;
      if (yy > 0) { double p = nf_series(yy, f->eps); if (p > 0 && R_FINITE(p)) ls = log(p) - 2.0*la; }
      for (int n = 0; n < Q; n++) {
        if (ls == LOGZERO) { f->lp[n] += LOGZERO; continue; }
        double v = (f->model == 1) ? (f->node[n] - f->vi[j]) : exp(f->node[n] - f->vi[j]);
        f->lp[n] += ls + v*a*xm - .5*v*v*t;
      }
    }
    s2 = exp(f->sd2A); Pval = par;
  } else {
    double Vp = (f->model == 1) ? par : exp(par);
    for (int j = 0; j < J; j++) {
      double yj = f->y[i + j*N];
      if (yj == -999) continue;
      double t = yj - f->xter[j], xm = f->x[i + j*N] - .5;
      double v = (f->model == 1) ? (Vp - f->vi[j]) : Vp/exp(f->vi[j]);
      for (int n = 0; n < Q; n++) {
        double ls = f->ls[j*Q + n];
        if (ls == LOGZERO) { f->lp[n] += LOGZERO; continue; }
        double a = f->a[j*Q + n];
        f->lp[n] += ls + v*a*xm - .5*v*v*t;
      }
    }
    s2 = exp(f->sd2V); Pval = par;   /* latent is v_p (D) or log theta_p (Q) */
  }
  double m = -INFINITY;
  for (int n = 0; n < Q; n++) { double z = f->lp[n] + f->logw[n]; if (z > m) m = z; }
  if (!(m > 0.5*LOGZERO) || !R_FINITE(m)) return -2.0*LOGFLOOR;   /* LLfacA/V: LL = 1e-300 */
  double s = 0.0;
  for (int n = 0; n < Q; n++) s += exp(f->lp[n] + f->logw[n] - m);
  /* normal prior on the latent variable: N(0, s2). (The original LLfacA/V used
     exp(-.5*P/sqrt(s2)), i.e. the deviation was neither squared nor divided by
     the variance.) */
  double lpost = m + log(s) - .5*log(2.0*M_PI*s2) - .5*Pval*Pval/s2;
  return -2.0*lpost;
}

/* Brent's minimizer on [lo, hi] (golden section + parabolic steps). */
static double brent_min(double lo, double hi, double (*fn)(double, void *), void *info, double tol)
{
  const double c = 0.3819660112501051, eps = 1.4901161193847656e-08;
  double a = lo, b = hi, v = a + c*(b - a), w = v, x = v, d = 0.0, e = 0.0;
  double fx = fn(x, info), fv = fx, fw = fx;
  for (int it = 0; it < 200; it++) {
    double xm = .5*(a + b), tol1 = eps*fabs(x) + tol/3.0, tol2 = 2.0*tol1;
    if (fabs(x - xm) <= tol2 - .5*(b - a)) break;
    double p = 0, q = 0, r = 0, u;
    if (fabs(e) > tol1) {
      r = (x - w)*(fx - fv); q = (x - v)*(fx - fw); p = (x - v)*q - (x - w)*r; q = 2.0*(q - r);
      if (q > 0) p = -p; else q = -q;
      r = e; e = d;
    }
    if (fabs(p) >= fabs(.5*q*r) || p <= q*(a - x) || p >= q*(b - x)) {
      e = (x < xm) ? b - x : a - x; d = c*e;
    } else {
      d = p/q; u = x + d;
      if (u - a < tol2 || b - u < tol2) d = (x < xm) ? tol1 : -tol1;
    }
    u = (fabs(d) >= tol1) ? x + d : (d > 0 ? x + tol1 : x - tol1);
    double fu = fn(u, info);
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return x;
}

void FACscore(double *y, double *x, int *N, int *nit, double *ai, double *vi,
              double *sd2A, double *sd2V, double *ter, int *which, int *model,
              double *A, double *W, int *nq, double *epsAlt, double *start,
              double *halfwidth, double *tol, int *adapt, double *muO, double *sO,
              double *est, double *hess, int *conv)
{
  int n = *N, J = *nit, Q = *nq;
  double *xter = (double *)R_alloc(J, sizeof(double));
  for (int j = 0; j < J; j++) xter[j] = exp(ter[j]);
  /* integrate over the *other* latent variable */
  double s2O = exp(*which == 1 ? *sd2V : *sd2A), sdOther = sqrt(s2O);
  int deg = !(s2O > 0);

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    facinfo f;
    double *node = (double *)malloc((size_t)Q*sizeof(double));
    double *logw = (double *)malloc((size_t)Q*sizeof(double));
    f.y = y; f.x = x; f.ai = ai; f.vi = vi; f.xter = xter; f.node = node; f.logw = logw;
    f.N = n; f.nit = J; f.nq = Q; f.model = *model; f.which = *which;
    f.sd2A = *sd2A; f.sd2V = *sd2V; f.eps = *epsAlt;
    f.lp = (double *)malloc((size_t)Q*sizeof(double));
    f.ls = (double *)malloc((size_t)J*Q*sizeof(double));
    f.a  = (double *)malloc((size_t)J*Q*sizeof(double));
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic, 8)
    #endif
    for (int i = 0; i < n; i++) {
      f.i = i;
      if (*adapt && !deg) {
        double s = (sO[i] > sdOther || !(sO[i] > 0)) ? sdOther : sO[i];
        for (int q = 0; q < Q; q++) {
          node[q] = muO[i] + s*A[q];
          logw[q] = log(W[q]) + log(s) - .5*log(s2O) - node[q]*node[q]/(2.0*s2O) + .5*A[q]*A[q];
        }
      } else {
        for (int q = 0; q < Q; q++) { node[q] = deg ? 0.0 : sdOther*A[q]; logw[q] = log(W[q]); }
      }
      if (f.which == 2) {            /* series depends only on the boundary node: cache it */
        for (int j = 0; j < J; j++) {
          double yj = y[i + j*n];
          for (int q = 0; q < Q; q++) {
            double la = node[q] - ai[j], a = exp(la), ls = LOGZERO;
            if (yj != -999) {
              double yy = (yj - xter[j])/(a*a);
              if (yy > 0) { double p = nf_series(yy, f.eps); if (p > 0 && R_FINITE(p)) ls = log(p) - 2.0*la; }
            }
            f.ls[j*Q + q] = ls; f.a[j*Q + q] = a;
          }
        }
      }
      double s0 = R_FINITE(start[i]) ? start[i] : 0.0;
      double lo = s0 - *halfwidth, hi = s0 + *halfwidth;
      double b = brent_min(lo, hi, fac_obj, &f, *tol);
      est[i] = b;
      /* not too small: the truncated series makes the objective step by ~1e-6 */
      /* where the number of terms changes, which a tiny h turns into garbage. */
      double h = 1e-2;
      hess[i] = (fac_obj(b + h, &f) - 2.0*fac_obj(b, &f) + fac_obj(b - h, &f))/(h*h);
      conv[i] = (b - lo < 1e-3*(*halfwidth) || hi - b < 1e-3*(*halfwidth)) ? 1 : 0;
    }
    free(f.lp); free(f.ls); free(f.a); free(node); free(logw);
  }
}
