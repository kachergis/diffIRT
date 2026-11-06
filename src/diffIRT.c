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
               * exp(-0.5 * Ap[i] / sqrt(exp(*sd2A)));
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
    else
      LL[i] = tmp2 / sqrt(2.0 * M_PI * exp(*sd2V))
               * exp(-0.5 * Vp[i] / sqrt(exp(*sd2V)));
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
  int xi = 0;
  double p0 = p / (aal*aal) * exp(exl*aal*(xi-.5) - .5*exl*exl*(*y-*ter));
  xi = 1;
  double p1 = p / (aal*aal) * exp(exl*aal*(xi-.5) - .5*exl*exl*(*y-*ter));
  *ou = p0 + p1;
}


/* ------------------------------------------------------------------ */
/* mdiffQQ  (unchanged)                                                 */
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
    D[i] = (tmp == 0.0 || !R_FINITE(tmp)) ? -999.0 : tmp;
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
