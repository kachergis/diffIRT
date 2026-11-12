# diffIRT 1.5.9000 (development version)

Changes since diffIRT 1.5 (CRAN, 2015). Model fits are 20–60x faster with the
original (non-adaptive) quadrature. The new default, adaptive quadrature, costs about
twice that but removes a bias in `omega[gamma]`. Factor scores are about 1,000x faster,
and the `RespFit()` test now works for long tests. Several long-standing bugs are fixed,
some of which change published results. Those changes are listed first.

## Changes that affect results

Results from diffIRT 1.5 and earlier will differ from this version in the ways
below. Where the old behaviour can be restored, the setting is given.

| Change | Effect on results | Restore old behaviour |
|---|---|---|
| Adaptive Gauss–Hermite quadrature is the default | `omega[gamma]` about 15% larger (the old estimate was ~13% too small); SEs of `a[i]` 20–25% larger (they were too small) | `control = list(adaptive = FALSE)` |
| Series precision `eps` default `.01` → `1e-6` | −2LL about 2 higher (more accurate); estimates change by < .01 | `control = list(eps = .01)` |
| Ter starting values fixed; restart of runaway parameters | Some fits that previously stopped early now reach the optimum (−2LL 50–680 lower in simulations) | `control = list(restart = Inf)` disables the restart; the starting-value fix cannot be disabled (supply `start` instead) |
| Standard errors of `omega[gamma]`, `omega[theta]` | About 2x smaller (they were too large) | — |
| Standard errors in models with equality or fixed constraints | All SEs change (each was scaled by the wrong parameter) | — |
| `vcov()` | Twice as large (it was half the covariance matrix) | — |
| `factest()` prior | Factor scores change, mostly for theta; scores of about 8% of persons no longer collapse to theta ≈ 0 in the Q model | — |
| `factest()` uses each person's adaptive grid | Scores and SEs more accurate: with 20 items, theta scores moved by up to 0.19 and SEs by up to 37% relative to a 61-point reference; now within ~1e-3 | `factest(control = list(adaptive = FALSE))` |
| `RespFit()` | Mr changes for the Q model; df changes for both models; p-values change | — |
| Input validation | Data that were silently accepted now give an error (responses other than 0/1/NA, non-positive RTs) | Recode the data |

### Details

* **Adaptive quadrature (`control$adaptive`, default `TRUE`).** With 20 items, each
  person's posterior for log boundary separation has an SD of about 0.1. That is
  narrower than the spacing of a 7-point Gauss–Hermite grid. As a result, the
  non-adaptive rule underestimated `omega[gamma]` by about 13% at the default `nq = 7`,
  regardless of N, so 95% intervals almost never covered the true value at N ≥ 500. It
  also understated the SEs of `a[i]` by 20–25%.

  Each person now gets a grid centred and scaled on their own posterior for boundary
  separation and drift. In a recovery simulation (20 items, D and Q models,
  N = 500 and 1000), bias is about 0 for every parameter. SE/SD ratios are 0.92–1.12 and
  coverage is .88–.98. At `nq = 7`, −2LL is within 0.004 of a 61-point grid, against 2–7
  units for the non-adaptive rule. It costs about twice the time of a non-adaptive fit.
  Use at least `nq = 5`.

* **Series precision (`control$eps`, default `1e-6`, was `.01`).** The truncated
  Navarro–Fuss series changes its number of terms with the RT and parameters, so with
  `eps = .01` −2LL jumped by up to about 1e-3. The optimizer occasionally stopped at
  such a jump and still reported convergence: in 16 simulated fits, one model Q fit
  per gradient type stopped 147–550 −2LL units short. With `1e-6`, all reached the
  optimum at about the same speed. `factest()` now defaults to the `eps` used to fit
  the model.

* **Fits that stopped early.**
  * **Ter starting values.** `get.start` replaced a negative EZ estimate of Ter by 1,
    which can exceed an item's fastest RT. Those responses then have zero density, and
    the first optimizer step threw that Ter to about 0, where the gradient on the log
    scale vanishes. optim then reported convergence. The only visible sign was a NaN
    standard error with `se = TRUE`. Ter now starts at half the item's fastest RT,
    including for items that share a constrained Ter.
  * **Restarts.** The optimizer could also drift until a variance or item drift ran off
    on the log scale (e.g. `omega[theta]` → 0) and stop there. `diffIRT()` now
    restarts, at most twice, when a free parameter ends more than `control$restart` log
    units (default 5) from its starting value. It resets those parameters and keeps the
    restart only if −2LL improves.
  * **Effect.** In 32 simulated fits (20 items, N = 500 and 1000), 9 improved by
    50–680 −2LL units and none got worse. **Published fits may have stopped early
    wherever an item's EZ Ter estimate was negative.** Refitting with this version and
    comparing −2LL is the easiest check.

* **Standard errors.**
  * `omega` is reported as exp(log(omega²)/2), so Var(omega) = (omega/2)²·Var(log
    omega²). The code used (omega/2)·Var, which made the SEs about 2x too large.
  * With equality or fixed constraints, the delta-method transformation indexed the full
    parameter vector, while the variances are indexed per free parameter. Each SE was
    therefore scaled by some other parameter. In the `brightness` example, drift SEs
    were 2.5x too small and the Ter and `omega[theta]` SEs 1.6–2.3x too large.
  * `vcov()` inverted the Hessian of −2LL directly. The covariance matrix is
    `solve(hessian / 2)`.

* **`factest()` prior.** The person prior in `LLfacA`/`LLfacV` was
  `exp(-.5 * P / sqrt(s2))`. A normal density needs `exp(-.5 * P^2 / s2)`. It is now a
  normal prior on log gamma and on theta (D model) or log theta (Q model). In
  simulation (N = 200, 20 items), theta's correlation with the true values rose from
  .854 to .870 (D) and from .554 to .614 (Q). Individual theta scores move by up to
  about 0.7 (D) and 0.9 (Q). Gamma scores are essentially unchanged.

* **`factest()` with adaptive grids.** When the model was fitted with adaptive
  quadrature (the default), `factest()` integrates over the other latent variable on
  each person's adaptive grid and starts each search at the person's posterior mean.
  Compared with a 61-point grid (500 persons, 20 items), the scores agree to about
  1e-4 and the SEs to 0.3%. With the non-adaptive 7-point grid, theta scores were off
  by up to 0.19 and SEs by up to 37%. The cost is unchanged.

* **`RespFit()` (Mr test).**
  * **Q-model derivatives.** `fitDeriv` tested `toupper(model) == "Q"`, which is never
    true for the numeric model code 2, so the item drifts stayed on the log scale and
    the Q-model derivatives were wrong.
  * **df formula.** The degrees of freedom subtracted 2 more than the number of
    parameters. In the D model, Mr behaved like χ²(18) but was compared with χ²(16)
    (8 items), so a correct model was rejected 10% of the time at α = .05.
  * **How df is now defined.** Derivatives are taken with respect to the free parameters,
    respecting constraints. Mr uses only the parameter combinations that responses
    identify: the drift scale is not identified in the D model, and only `a[i] v[i]`
    and omega[gamma]² + omega[theta]² are in the Q model. df is the number of moments
    minus that count.
  * **Calibration.** With data simulated from the model (8 items, N = 500, 100
    replications), 6% of fits are rejected at .05 for both models, and the mean of Mr
    equals its df.

## Performance

* **Likelihood and gradient (`LLdiff`, `DERdiff`).**
  * **Fewer series evaluations.** The Navarro–Fuss series depends only on the boundary
    node, so it is now evaluated `nq` instead of `nq²` times per response.
  * **Log space.** The drift part is computed in log space, and the quadrature sum uses
    log-sum-exp, so long tests no longer underflow.
  * **OpenMP.** The code runs in parallel over subjects, with a deterministic sum, so
    results do not depend on the number of threads.
  * **Timings** (N = 1000, 20 items, `nq = 7`): likelihood 47 ms → 9.2 ms (1 thread) →
    1.4 ms (14 threads). The old forward-difference gradient went from 2.9 s to 62 ms
    (1 thread) and 8.8 ms (14 threads).
  * **Full fit** (N = 200, 20 items): 63 s → 2.9 s (1 thread) → 1.1 s (14 threads).
* **Analytic gradient (`control$grad`, default `"analytic"`).** It costs about one
  likelihood evaluation, against about six for the cached forward-difference gradient.
  It matches central differences to about 1e-8 relative error, against about 1e-5
  before. Full fits are about 2–2.5x faster again: N = 1000, 20 items, 8 datasets per
  model, D 18.0 → 7.7 s and Q 23.4 → 11.6 s in total. `grad = "numeric"` keeps the
  forward difference (non-adaptive fits only).
* **`factest()`.** The scores are separable across persons, so each person is now scored
  with a one-dimensional Brent search in C, in parallel. The old method was a joint
  BFGS over all N scores with a numerical gradient, whose cost grew as N². N = 200: D
  8.2 s → 0.009 s, Q 16.8 s → 0.004 s. N = 1000 with `se = TRUE`: 0.015 s.
* **`RespFit()`.** It no longer enumerates all 2^`nit` response patterns. Everything is
  computed from joint probabilities of at most 2 × `order` items. 14 items: 44 s →
  0.4 s. 20 items: 1.3 s. 40 items: 11 s. It was previously infeasible beyond about 14
  items.
* **`QQdiff()`.** The predicted RT distribution of an item is integrated once on a grid
  (500 panels of 8-point Gauss–Legendre, narrow near Ter), and all quantiles come from
  inverting its cubic Hermite interpolant. Previously each quantile called `uniroot()` on
  `integrate()`. N = 1000, 5 items, one quantile per person: about 30 s → 0.2 s. The
  quantiles agree with a high-precision reference to about 1e-8; the old ones were
  accurate only to about 1e-4 (the default `uniroot()` tolerance). The default number of
  quantiles is now `min(N, 100)` instead of N
  ([#12](https://github.com/kachergis/diffIRT/issues/12)).
* **The `nq` ramp-up** (e.g. `nq = c(1, 3, 5, 7)`) no longer pays off. The optimizer
  restarts at each stage, so with the faster likelihood a full ramp is slower than
  `nq = 7` alone. `nq = c(1, 7)` can be slightly faster.

## Bug fixes

* **Apple Silicon hang.** When the optimizer tried an extreme step, (y − Ter)/a² could be
  0 or NaN. `(int)ceil(Inf)` then saturates to `INT_MAX` on ARM64, so the series loop
  ran about 2×10⁹ times (about 16 s per density), and fits could take hours. On x86 the
  same cast gives `INT_MIN` and the loop is skipped, which is probably why this went
  unnoticed. Fixed in `ddiff`, `ddiffFacV`, `ddiffQQ` and the new code.
* **`factest(se = TRUE)`** always failed with an error. It now works, and the gamma SE
  column no longer copies the theta SEs.
* **Input validation in `diffIRT()`.** The following now stop with a clear error before
  any computation: responses other than 0, 1 or NA; non-positive RTs; and persons with no
  observed responses. The old check for all-missing persons
  (`length(complete.cases(rt)) != N`) could never trigger, so such data crashed with an
  uninformative error.
* **`anova()`** compared `object$x`, which does not exist, so its "different dataset"
  warning never fired. It now compares the responses and RTs.
* **`RespFit()` with missing data.** Missing responses were counted as incorrect in the
  observed moments, which biased them downward. Each moment is now computed from the
  persons who answered all its items, and the covariance of two moments is scaled by
  the number who answered the items of both (valid under MCAR). Results without
  missing data are unchanged ([#11](https://github.com/kachergis/diffIRT/issues/11)).
  With 10% and 20% MCAR missingness (8 items, N = 500, 100 replications;
  `sim/respfit_missing.R`), mean Mr is 18.2–18.6 for df 19 (D) and 26.1–27.0 for
  df 27 (Q), with 2–8% rejected at .05. Counting missing as incorrect rejected 31–75%.
  `print` and `summary` show the order actually used instead of always "Order 2".
* **Predicted RT density (`QQdiff()`).** The density computed exp(v·a/2) before
  multiplying by exp(−v²t/2), which overflows to Inf on the extreme nodes of a fine
  grid. It returned −999 where the density is 0, which `pdiffIRT` would then integrate
  as a value. The density is now computed on the log scale and returns 0. With the
  default 7-point grid neither problem showed up in practice: densities integrate to 1,
  as before. The histogram's predicted density was drawn from log Ter instead of Ter.
* **`benchmark.R`** passed `-fopenmp` only to the linker, so its "OpenMP" build ran
  serially.

## New options

* `diffIRT(control = list(...))`:
  * `adaptive`: adaptive quadrature (default `TRUE`).
  * `adapttol`, `maxadapt`: re-centring tolerance in −2LL (default `.01`) and maximum
    number of re-centrings (default 10).
  * `grad`: `"analytic"` (default) or `"numeric"`.
  * `restart`: restart threshold in log units (default 5; `Inf` disables).
  * `threads`: number of OpenMP threads (default `NULL`, the OpenMP default). Use
    `threads = 1` inside an outer parallel loop. Estimates do not depend on it
    ([#6](https://github.com/kachergis/diffIRT/issues/6)).
* `factest(control = list(...))`:
  * `tol`: tolerance of the per-person search (default `1e-8`).
  * `halfwidth`: half-width of the search interval (default `max(3, 10 × SD)`). A
    warning is given when an estimate ends at the edge of the interval.
  * `threads`: number of OpenMP threads (default: the value used for the fit).
* **EAP scores and plausible values** in `factest()`
  ([#8](https://github.com/kachergis/diffIRT/issues/8)):
  * `type = "EAP"`: posterior means, and posterior SDs with `se = TRUE`, on a 21 × 21
    adaptive Gauss–Hermite grid per person (`control$nq.eap`). They agree with a
    brute-force posterior on a 241 × 241 grid to about 1e-7, and take about 0.03 s for
    400 persons.
  * `type = "PV"`, `npv = 5`: plausible values, as a list of `npv` score matrices. Each
    value is drawn by importance resampling (`control$pv.draws = 200` candidates from a
    t(5) proposal centred on the person's posterior). When the fit has standard errors,
    each set first draws the item and population parameters from their asymptotic
    normal distribution (`control$pv.items`), redrawing any set in which a Ter exceeds
    the item's fastest RT.
  * In simulation (`sim/pv_variance.R`: 100 replications per model, N = 300, 10 items,
    10% missing), the variance of the plausible values across persons matches the
    population variance: for example 0.249 against 0.25 for theta in the D model, where
    MAP and EAP scores give 0.138. Q-model theta: plausible values 0.095 against 0.09;
    MAP 0.017.
  * `type = "MAP"` (the default) is unchanged.
  * New C entry point `LLpoints`: each person's log-likelihood at arbitrary latent points.
* **Fit statistics for response times and persons**
  ([#10](https://github.com/kachergis/diffIRT/issues/10)):
  * `RTFit(object, bins = 10)` gives a chi-square test per item and over all items. It
    compares observed and expected proportions of response (0/1) × RT bin, with bin
    edges at the item's predicted RT deciles. The residual covariance accounts for
    parameter estimation through Cov(p, θ̂) = ΔV for the maximum likelihood estimate,
    so no bootstrap or refitting is needed. df = 2 × bins − 1 per item.
  * `PersonFit(object, nsim = 1000)` gives each person's marginal log-likelihood with a
    Monte Carlo p-value. The reference is persons simulated from the fitted model with
    the same items answered. It also reports a standardized `lz`.
  * Calibration and power (`sim/fit_statistics.R`: 100 replications, N = 500, 10 items,
    10% missing; rejection rates at .05):

    | Condition | `RTFit` item | `RTFit` overall | `PersonFit` |
    |---|---|---|---|
    | Null, D / Q | .052 / .052 | .06 / .09 | .049 / .048 |
    | Item 1 RT +0.3 s for 30% of persons, D / Q | .58 / .52 (item 1) | .25 / .33 | — |
    | 5% inattentive on 5 items (RT ×2, random responses), D / Q | — | .42 / .47 | .64 / .63 |
    | 5% rapid guessing on 5 items, D / Q | — | .49 / .57 | .10 / .09 |

    Rapid guessing is hard to detect at the person level. The guessers pull the
    estimated Ter down, which masks them: with the true item parameters, 88% were
    flagged. Fast responses at chance accuracy on *every* item fit the model as a person
    with a small boundary, so they show up as a low `gamma[p]` rather than as misfit. A
    purification step (refit without the flagged persons) raised the Type I error to
    9% and power only to 24%, so it was not added.
  * New C entry point `ddiffBase`: node-wise first-passage density.
* `simdiff()`:
  * `omega`: standard deviations of the person parameters (log gamma, and theta or log
    theta), when they are generated.
  * `rho`: correlation between the two person latent variables (note that `diffIRT()`
    itself still assumes independence; see [#2](https://github.com/kachergis/diffIRT/issues/2)).
  * `miss`: proportion of responses set to `NA` completely at random.

  With the default arguments the random numbers are drawn exactly as before, so
  existing seeds reproduce the same data.
* New methods for fitted models
  ([#9](https://github.com/kachergis/diffIRT/issues/9)):
  * `predict()`: expected P(correct) and mean RT per item, for the population or for
    given person parameters.
  * `fitted()`: the same at each person's factor scores.
  * `simulate()`: new datasets from the fitted model, optionally with the same missing
    pattern, for parametric bootstraps and posterior predictive checks.
* The fit object stores the adaptive grids (`$aq`). The per-person log-likelihoods
  (`$subjLL`) use them.

## Testing and validation

* **New `testthat` suite** (`tests/testthat`, about 6 s; the slower tests skip on CRAN):
  * per-person likelihoods match the original C code, using stored fixtures;
  * the ARM64 guard returns quickly;
  * results are identical with 1 and 4 OpenMP threads;
  * analytic and adaptive gradients match central differences, with missing data,
    constraints, fixed parameters and capped grids;
  * adaptive-quadrature special cases;
  * parameter recovery;
  * `vcov()` and delta-method SEs on the constrained `brightness` model;
  * Ter starting values and the rescue of a collapsed fit;
  * `factest()` scores are stationary points;
  * `RespFit()` moments match full pattern enumeration, and its df and calibration are
    right;
  * input validation and `anova()`.

  Reintroducing the old omega-SE formula or a wrong gradient term makes the suite fail.
  `R CMD check`: Status OK.
* **Parameter recovery simulation:** `sim/recovery.R` and `sim/summarize.R` (bias, RMSE,
  SE/SD ratio, coverage, factor-score recovery). They are excluded from the built
  package.

## Internal changes

* **New C routines:** `DERdiffAn` (analytic gradient); `LLdiffAQ`, `DERdiffAQ` and
  `AQmoments` (adaptive quadrature); `FACscore` (factor scores).
* **New R function:** `aqMoments`.
* **Removed:** `fitPred`, `fitDeriv` and `fitOrders`, which only the old `RespFit()`
  used.
* **Packaging:** `testthat` is in Suggests; `.Rbuildignore` excludes development files
  (`src_orig/`, `benchmark.R`, `sim/`, `extra/`, PDFs).
* **OpenMP setup.** `src/Makevars` uses `$(SHLIB_OPENMP_CFLAGS)`. On macOS, OpenMP
  needs a compiler that supports it (e.g. Homebrew LLVM via `~/.R/Makevars`); with
  Apple clang the package builds and runs serially. `omp.h` is no longer included.

## Known limitations and planned work

Tracked as issues on GitHub:

* `QQdiff()` and `RespFit()` integrate population-level quantities over the prior,
  where the non-adaptive 7-point grid is accurate: the predicted RT density is within
  0.1% of a 61-point grid, and the `RespFit()` moments within 0.001
  ([#7](https://github.com/kachergis/diffIRT/issues/7), resolved together with
  adaptive grids in `factest()`).
* The person boundary (gamma) and drift (theta) are assumed independent. A correlation
  between them is not yet supported ([#2](https://github.com/kachergis/diffIRT/issues/2)).
* `simdiff()`'s rejection sampler can fail (stop at `max.iter`) for extreme
  person–item combinations where |v·a| is very large (responses almost certainly
  correct). This is rare with the default parameter distributions. The sampler is
  fast enough in R (0.2 s for 1000 persons × 20 items), so it was not ported to C
  ([#9](https://github.com/kachergis/diffIRT/issues/9)).
* Planned model extensions:
  * start-point bias ([#3](https://github.com/kachergis/diffIRT/issues/3));
  * between-trial variability in drift, start point and non-decision time
    ([#4](https://github.com/kachergis/diffIRT/issues/4));
  * person covariates and multiple groups ([#5](https://github.com/kachergis/diffIRT/issues/5)).

`src_orig/` holds the unmodified original C code, used for benchmarking and for the test
fixtures. It is excluded from the built package.
