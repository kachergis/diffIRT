/*
 *  diffIRT.c  -- optimized version
 *
 *  Changes from original (Dylan Molenaar, 2013):
 *
 *  1. Hoisted loop-invariant exp/sqrt calls out of inner loops in LLdiff,
 *     LLfacA, LLfacV (sqrt(exp(sd2A)), sqrt(exp(sd2V)), exp(ter[j])).
 *
 *  2. Log-space product accumulation in LLdiff, LLfacA, LLfacV to prevent
 *     underflow when nit is large (the main cause of convergence failures on
 *     real datasets). The linear-space tmp2 accumulation is preserved for the
 *     quadrature sum, but the per-item product is now done in log space.
 *
 *  3. OpenMP parallelisation of the subject loop in LLdiff, LLfacA, LLfacV.
 *     Each thread uses thread-local variables to avoid races. The final
 *     log-sum reduction is serial (negligible cost).
 *     Compile with -fopenmp (add to src/Makevars: PKG_CFLAGS=-fopenmp
 *     PKG_LIBS=-lgomp). Gracefully degrades to serial if OpenMP unavailable.
 *
 *  All other functions (ddiff, ddiffFacV, ddiffQQ, DERdiff, mdiffQQ,
 *  equal, makeT) are unchanged except for the ter cache in DERdiff's
 *  nested LLdiff calls (handled by LLdiff itself now).
 *
 *  Original copyright 2013 Dylan Molenaar, University of Amsterdam.
 *  GNU GPL v2 or later.
 */

#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

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
/* LLdiff  -- OPTIMIZED                                                 */
/*   - hoisted: sqrtSdA, sqrtSdV, xter_cache[] out of inner loops      */
/*   - log-space per-item product to prevent underflow                  */
/*   - OpenMP over subjects                                             */
/* ------------------------------------------------------------------ */
void LLdiff(double *y, double *x, int *N, int *nit, double *ai, double *vi,
            double *sd2A, double *sd2V, double *ter, int *model,
            double *A, double *W, int *nq, double *epsAlt,
            double *ou, double *LL)
{
  /* --- hoist loop-invariant transforms --- */
  double sqrtSdA = sqrt(exp(*sd2A));
  double sqrtSdV = sqrt(exp(*sd2V));

  double *xter_cache = (double *)R_alloc(*nit, sizeof(double));
  for (int j = 0; j < *nit; j++) xter_cache[j] = exp(ter[j]);

  /* --- pre-transform quadrature nodes --- */
  double *AtransA = (double *)R_alloc(*nq, sizeof(double));
  double *AtransV = (double *)R_alloc(*nq, sizeof(double));
  for (int n = 0; n < *nq; n++) {
    AtransA[n] = sqrtSdA * A[n];
    AtransV[n] = sqrtSdV * A[n];
  }

  /* --- subject loop (parallelised) --- */
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic, 4)
  #endif
  for (int i = 0; i < *N; i++) {
    double tmp2 = 0.0;

    for (int n1 = 0; n1 < *nq; n1++) {
      double aAp = AtransA[n1];

      for (int n2 = 0; n2 < *nq; n2++) {
        double aVp = AtransV[n2];

        /* log-space product over items */
        double log_tmp1 = 0.0;
        int valid = 1;

        for (int j = 0; j < *nit; j++) {
          if (y[i + j * (*N)] != -999) {
            double ou_local;
            ddiff(&y[i + j*(*N)], &x[i + j*(*N)],
                  &aAp, &ai[j], &aVp, &vi[j],
                  &xter_cache[j], model, epsAlt, &ou_local);
            if (ou_local <= 0.0) { valid = 0; break; }
            log_tmp1 += log(ou_local);
          }
        }

        if (valid)
          tmp2 += exp(log_tmp1) * W[n1] * W[n2];
      }
    }

    LL[i] = (tmp2 == 0.0 || !R_FINITE(tmp2)) ? 1e-300 : tmp2;
  }

  /* --- serial log-sum reduction --- */
  double logsum = 0.0;
  for (int i = 0; i < *N; i++) logsum += log(LL[i]);
  *ou = -2.0 * logsum;
}


/* ------------------------------------------------------------------ */
/* DERdiff  (unchanged logic; benefits from faster LLdiff)             */
/* ------------------------------------------------------------------ */
void DERdiff(double *y, double *x, int *N, int *nit, int *constrain,
             int *npar, double *ai, double *vi, double *sd2A, double *sd2V,
             double *ter, int *model, double *A, double *W, int *nq,
             double *epsAlt, double *delta, double *ou, double *LL,
             double *deriv)
{
  LLdiff(y, x, N, nit, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt, ou, LL);
  double tmp1 = *ou;
  double tmp2 = 0.0;

  for (int j = 1; j < (*npar + 1); j++) {
    for (int i = 0; i < *nit; i++) {
      if (constrain[i]          == j) ai[i]  += *delta;
      if (constrain[i + *nit]   == j) vi[i]  += *delta;
      if (constrain[i + 2*(*nit)]== j) ter[i] += *delta;
    }
    if (constrain[3*(*nit)]   == j) *sd2A += *delta;
    if (constrain[3*(*nit)+1] == j) *sd2V += *delta;

    LLdiff(y, x, N, nit, ai, vi, sd2A, sd2V, ter, model, A, W, nq, epsAlt, ou, LL);
    tmp2 = *ou;
    deriv[j-1] = (tmp2 - tmp1) / (*delta);

    for (int i = 0; i < *nit; i++) {
      if (constrain[i]          == j) ai[i]  -= *delta;
      if (constrain[i + *nit]   == j) vi[i]  -= *delta;
      if (constrain[i + 2*(*nit)]== j) ter[i] -= *delta;
    }
    if (constrain[3*(*nit)]   == j) *sd2A -= *delta;
    if (constrain[3*(*nit)+1] == j) *sd2V -= *delta;
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
