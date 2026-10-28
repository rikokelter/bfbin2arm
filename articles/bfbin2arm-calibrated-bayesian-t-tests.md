# Calibrated Bayesian One-Sample and Paired t-Tests for Continuous Endpoints in Clinical Trials

``` r

library(bfbin2arm)

# Use default controls for reproducible vignette results.
# For interactive exploratory work, users may supply a faster control list.
control_fast <- list(
  n_search = "two_stage",
  n_step_stage1 = 10L,
  rel_tol_m1 = 1e-6,
  subdivisions_m1 = 100L,
  grid_n = 121L,
  t_max = 6,
  t_max_cap = 16,
  root_tol = 1e-8,
  verify = FALSE
)

control_full_fast <- list(
  n_search = "full_fast",
  rel_tol_m1 = 1e-6,
  subdivisions_m1 = 100L,
  grid_n = 121L,
  t_max = 6,
  root_tol = 1e-8,
  verify = FALSE,
  progress = TRUE
)

vignette_control <- list(
  n_search = "full",
  rel_tol_m1 = 1e-9,
  subdivisions_m1 = 800L,
  grid_n = 801L,
  t_max = 12,
  t_max_cap = 40,
  root_tol = 1e-12,
  verify = TRUE,
  verify_rel_tol_m1 = 1e-10,
  verify_subdivisions_m1 = 1200L,
  verify_grid_n = 1201L,
  verify_t_max = 15,
  verify_root_tol = 1e-12,
  progress = TRUE
)
```

## Overview

This vignette illustrates deterministic calibration of Bayesian
one-sample and paired $`t`$-test designs for continuous clinical
endpoints. The examples are prospective, trial-inspired design
illustrations; they do not reproduce the original trials’ analyses or
numerical results.

The function
[`calibrate_bf_ttest()`](https://rikokelter.github.io/bfbin2arm/reference/calibrate_bf_ttest.md)
searches the requested range of sample sizes and returns the smallest
feasible integer sample size. It evaluates evidence using the Bayes
factor $`BF_{01}`$, where smaller values favour $`H_1`$ and larger
values favour $`H_0`$.

The package supports three hypotheses:

``` math
H_0:\delta=0 \quad \text{versus} \quad H_1:\delta\ne0,
```

``` math
H_0:\delta=0 \quad \text{versus} \quad H_1:\delta>0,
```

and

``` math
H_0:\delta\le0 \quad \text{versus} \quad H_1:\delta>0.
```

Here, $`\delta`$ is the standardized mean effect. For directional
examples, the endpoint is coded so that positive values represent
clinical benefit.

## Example 1: Two-sided point-null design

This example is inspired by a phase-2 COPD study involving ivacaftor,
compare Vijaykumar et al. (2025). We consider a hypothetical single-arm
proof-of-concept study whose endpoint is change in prebronchodilator
forced expiratory volume in 1 second (FEV1) percent predicted at Day 84.

Suppose that a 5-percentage-point improvement is clinically relevant and
that the planning standard deviation is 12.5 percentage points. The
corresponding standardized effect is

``` math
\delta_A = \frac{5}{12.5} = 0.40.
```

We use the two-sided point-null hypothesis

``` math
H_0:\delta=0
\quad\text{versus}\quad
H_1:\delta\ne0.
```

The analysis prior under $`H_1`$ is $`N(0, 0.5^2)`$. This prior defines
the Bayes factor used at analysis; it is not the same as the design
prior used to evaluate operating characteristics.

The design prior under $`H_0`$ is a point mass at $`\delta=0`$. Under
$`H_1`$, we use a truncated normal design prior centred at the
anticipated beneficial effect, $`TN(0.40, 0.20^2; [0,\infty))`$. This
represents uncertainty about the size of a beneficial effect while
excluding harmful effects from the beneficial-effect planning scenario.

``` r

pulmonary_normal_dp <- calibrate_bf_ttest(
  design = "one.sample",
  hypothesis = "two.sided.pointnull",

  # Search range for the number of evaluable participants
  n_min = 20,
  n_max = 200,

  # Evidence thresholds:
  # BF01 <= 0.10 is compelling evidence for H1;
  # BF01 >= 3 is compelling evidence for H0.
  k1 = 1 / 10,
  k0 = 3,

  # Analysis prior under H1: N(0, 0.5^2)
  analysis_prior_sd = 0.5,

  # Error probability for incorrectly obtaining evidence for H1
  # when the true standardized effect is exactly zero.
  freq_t1e = 0.05,
  bayes_t1e = 0.05,

  # Frequentist power at the clinically relevant point alternative.
  freq_power = list(
    main = list(effect = 0.40, target = 0.80)
  ),

  # Average probability of evidence for H1 under the H1 design prior.
  bayes_power = 0.80,

  # Probability of compelling evidence for H0 under the H0 design prior.
  # Including this argument also computes and stores the H0 evidence region.
  pce_h0 = 0.60,

  # Design priors used only for operating-characteristic calculations.
  design_prior_h0 = list(type = "point", effect = 0),
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = control_full_fast
)

print(pulmonary_normal_dp)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: two.sided.pointnull 
    n*: 161 (df = 160 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: N(0, 0.5^2) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, Inf]
    Evidence for H1 region(s) on t-scale:
            lower     upper
    [1,]     -Inf -2.951809
    [2,] 2.951809       Inf
    Evidence for H0 region(s) on t-scale:
             lower    upper
    [1,] -1.248312 1.248312

    Constraints:
          constraint        type
            freq_t1e frequentist
           bayes_t1e    Bayesian
         bayes_power    Bayesian
              pce_h0    Bayesian
     freq_power_main frequentist
                    estimand    estimate target
       P_delta=0(BF01 <= k1) 0.003634714   0.05
     P_design,H0(BF01 <= k1) 0.003634714   0.05
     P_design,H1(BF01 <= k1) 0.800167910   0.80
     P_design,H0(BF01 >= k0) 0.786259448   0.60
     P_delta=0.4(BF01 <= k1) 0.982124193   0.80
     direction      slack satisfied
            <= 0.04636529      TRUE
            <= 0.04636529      TRUE
            >= 0.00016791      TRUE
            >= 0.18625945      TRUE
            >= 0.18212419      TRUE

The printed result reports the selected sample size $`n^*`$, the
decision regions on the observed-$`t`$ scale, and the attained values of
every active constraint. For this configuration, the calibrated result
is $`n^*=161`$ ($`df=160`$); the evidence-for-$`H_1`$ region is
approximately $`t\le -2.95`$ or $`t\ge 2.95`$.

``` r

plot(pulmonary_normal_dp)
```

![Figure 1: Output of the plot function for a one-sample Bayesian
two-sided point-null t-test design based on full calibration. Both
Bayesian and frequentist power and type-I-error are calibrated
jointly.](figures/bfbin2arm-onesample-t-test-ex1.png)

Figure 1: Output of the plot function for a one-sample Bayesian
two-sided point-null t-test design based on full calibration. Both
Bayesian and frequentist power and type-I-error are calibrated jointly.

The frequentist type-I error is the probability, at $`\delta=0`$, of
obtaining $`BF_{01}\le 0.10`$, that is, compelling evidence for $`H_1`$
when the point null is true. The Bayesian type-I error averages that
same probability over the $`H_0`$ design prior. Because the $`H_0`$
design prior here is a point mass at zero, these two values coincide.

The operating-characteristic plot displays their shared curve as
“Frequentist type-I error” rather than drawing two coincident curves and
repeating the legend entry.

The frequentist power constraint evaluates
$`P_{\delta=0.40}(BF_{01}\le 0.10)`$. In contrast, `bayes_power`
averages $`P_{\delta}(BF_{01}\le 0.10)`$ over the continuous $`H_1`$
design prior. Thus, the two power requirements answer different planning
questions.

The `pce_h0` constraint evaluates the probability of obtaining
compelling evidence for the null, $`BF_{01}\ge3`$, under the $`H_0`$
design prior. This implies that a minimum probability of obtaining a
Bayes factor in favour of $`H_0`$ of at least size $`3`$ or larger is
fixed, when $`H_0`$ is true. In this example, it also causes the
function to calculate the central evidence-for-$`H_0`$ region. The
interval is shown by the printed object and by the blue-shaded region in
the decision-function plot in the top left.

The decision panel displays $`\log(BF_{01})`$. The orange dashed line
denotes $`\log(k_1)`$, so values below it support $`H_1`$; the blue
dashed line denotes $`\log(k_0)`$, so values above it support $`H_0`$.
The shaded regions represent the stored decision regions: orange for
evidence for $`H_1`$, blue for evidence for $`H_0`$, and grey for
inconclusive results. The sample-size panel shows the operating
characteristics throughout the search and marks the selected $`n^*`$.

### A brief comparison of calibration modes

The first example above calibrated both frequentist and Bayesian
operating characteristics such as power and type-I-error jointly. Also,
the probability of compelling evidence for $`H_0`$ was calibrated. Next,
we compare the different calibration modes available. We can pick
whatever combination we like as the calibration target, but the most
prominent choices might include:

- Frequentist calibration
- Bayesian calibration
- Hybrid calibration
- Full calibration

We use identical arguments including our design and analysis priors for
all calibrations below, as specified in the following code chunk:

``` r

library(bfbin2arm)

vignette_control <- resolve_control_bf_ttest()

common_args <- list(
  design = "one.sample",
  hypothesis = "two.sided.pointnull",

  n_min = 20,
  n_max = 200,

  # BF01 <= 0.10: evidence for H1
  # BF01 >= 3: available as a reference threshold, but no PCE constraint here
  k1 = 1 / 10,
  k0 = 3,

  # Analysis prior under H1: N(0, 0.5^2)
  analysis_prior_sd = 0.5,

  # H0 design prior: exact null
  design_prior_h0 = list(type = "point", effect = 0),

  # H1 design prior: uncertain beneficial effect
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = control_full_fast
)
```

#### Frequentist calibration

This design requires frequentist type-I error at most 0.05 and
frequentist power at $`\delta=0.40`$ of at least 0.80.

``` r

pulmonary_freq <- do.call(
  calibrate_bf_ttest,
  c(
    common_args,
    list(
      freq_t1e = 0.05,
      freq_power = list(
        main = list(effect = 0.40, target = 0.80)
      )
    )
  )
)

print(pulmonary_freq)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: two.sided.pointnull 
    n*: 88 (df = 87 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: N(0, 0.5^2) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, Inf]
    Evidence for H1 region(s) on t-scale:
            lower     upper
    [1,]     -Inf -2.897943
    [2,] 2.897943       Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
          constraint        type                estimand
            freq_t1e frequentist   P_delta=0(BF01 <= k1)
     freq_power_main frequentist P_delta=0.4(BF01 <= k1)
        estimate target direction        slack satisfied
     0.004751998   0.05        <= 0.0452480024      TRUE
     0.800302685   0.80        >= 0.0003026846      TRUE

#### Bayesian calibration

This design requires Bayesian type-I error at most 0.05 and Bayesian
power of at least 0.80 under the specified $`H_1`$ design prior.

``` r

pulmonary_bayes <- do.call(
  calibrate_bf_ttest,
  c(
    common_args,
    list(
      bayes_t1e = 0.05,
      bayes_power = 0.80
    )
  )
)

print(pulmonary_bayes)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: two.sided.pointnull 
    n*: 161 (df = 160 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: N(0, 0.5^2) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, Inf]
    Evidence for H1 region(s) on t-scale:
            lower     upper
    [1,]     -Inf -2.951809
    [2,] 2.951809       Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
      constraint     type                estimand    estimate
       bayes_t1e Bayesian P_design,H0(BF01 <= k1) 0.003634714
     bayes_power Bayesian P_design,H1(BF01 <= k1) 0.800167910
     target direction      slack satisfied
       0.05        <= 0.04636529      TRUE
       0.80        >= 0.00016791      TRUE

#### Hybrid calibration

This hybrid design controls type-I error at the exact null and requires
average Bayesian power under the H1 design prior.

``` r

pulmonary_bayes_power_freq_t1e <- do.call(
  calibrate_bf_ttest,
  c(
    common_args,
    list(
      freq_t1e = 0.05,
      bayes_power = 0.80
    )
  )
)

print(pulmonary_bayes_power_freq_t1e)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: two.sided.pointnull 
    n*: 161 (df = 160 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: N(0, 0.5^2) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, Inf]
    Evidence for H1 region(s) on t-scale:
            lower     upper
    [1,]     -Inf -2.951809
    [2,] 2.951809       Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
      constraint        type                estimand    estimate
        freq_t1e frequentist   P_delta=0(BF01 <= k1) 0.003634714
     bayes_power    Bayesian P_design,H1(BF01 <= k1) 0.800167910
     target direction      slack satisfied
       0.05        <= 0.04636529      TRUE
       0.80        >= 0.00016791      TRUE

#### Full calibration

This activates all four error/power constraints, but deliberately does
not require a minimum probability of compelling evidence for $`H_0`$.

``` r

pulmonary_full_no_pce <- do.call(
  calibrate_bf_ttest,
  c(
    common_args,
    list(
      freq_t1e = 0.05,
      bayes_t1e = 0.05,
      freq_power = list(
        main = list(effect = 0.40, target = 0.80)
      ),
      bayes_power = 0.80

      # Deliberately omit:
      # pce_h0 = ...
    )
  )
)

print(pulmonary_full_no_pce)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: two.sided.pointnull 
    n*: 161 (df = 160 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: N(0, 0.5^2) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, Inf]
    Evidence for H1 region(s) on t-scale:
            lower     upper
    [1,]     -Inf -2.951809
    [2,] 2.951809       Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
          constraint        type                estimand
            freq_t1e frequentist   P_delta=0(BF01 <= k1)
           bayes_t1e    Bayesian P_design,H0(BF01 <= k1)
         bayes_power    Bayesian P_design,H1(BF01 <= k1)
     freq_power_main frequentist P_delta=0.4(BF01 <= k1)
        estimate target direction      slack satisfied
     0.003634714   0.05        <= 0.04636529      TRUE
     0.003634714   0.05        <= 0.04636529      TRUE
     0.800167910   0.80        >= 0.00016791      TRUE
     0.982124193   0.80        >= 0.18212419      TRUE

The frequentist-only calibration selects $`n^*=88`$, because it requires
80% probability of evidence for $`H_1`$ only at the single planning
effect $`\delta=0.40`$. The Bayesian-power requirement averages this
probability over the H1 design prior, which includes smaller positive
effects; this produces the larger sample size $`n^*=161`$.

The Bayesian-only, Bayesian-power-plus-frequentist-type-I-error, and
full no-PCE calibrations select the same sample size in this
illustration. This occurs because the $`H_0`$ design prior is a point
mass at $`\delta=0`$, making the Bayesian and frequentist type-I error
quantities numerically identical, and because the attained frequentist
power at $`\delta=0.40`$ already exceeds its 0.80 target at $`n=161`$.
No calibration includes a PCE requirement. Hence, none requires a
minimum probability of compelling evidence for $`H_0`$.

## Example 2: Directional point-null design

Consider the same pulmonary endpoint when the confirmatory scientific
question is restricted to benefit. The directional point-null test
compares $`H_0:\delta=0`$ with $`H_1:\delta>0`$.

The analysis prior under $`H_1`$ is a normal distribution truncated to
positive standardized effects, $`N(0, 0.5^2)\mathbb{I}(\delta>0)`$.
Consequently, an adverse observed change does not provide evidence for
the beneficial alternative. The design prior under $`H_1`$ represents
uncertainty about the anticipated beneficial effect.

``` r

directional_common_args <- list(
  design = "one.sample",
  hypothesis = "greater.pointnull",

  # Search range for the number of evaluable participants.
  n_min = 20,
  n_max = 200,

  # BF01 <= 0.10: compelling evidence for benefit (H1).
  # BF01 >= 3: available as a decision threshold, but no PCE constraint is
  # imposed in the calibrations below.
  k1 = 1 / 10,
  k0 = 3,

  # Directional H1 analysis prior:
  # N(0, 0.5^2), truncated to delta > 0.
  analysis_prior_sd = 0.5,

  # H0 design prior: exact null effect.
  design_prior_h0 = list(type = "point", effect = 0),

  # H1 design prior: uncertain beneficial effect.
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = control_full_fast
)
```

### Frequentist calibration

This design controls frequentist type-I error at the exact null and
requires frequentist power of at least 80% at $`\delta=0.40`$.

``` r

pulmonary_directional_freq <- do.call(
  calibrate_bf_ttest,
  c(
    directional_common_args,
    list(
      freq_t1e = 0.05,
      freq_power = list(
        main = list(effect = 0.40, target = 0.80)
      )
    )
  )
)

print(pulmonary_directional_freq)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: greater.pointnull 
    n*: 75 (df = 74 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 2.601071   Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
          constraint        type                estimand
            freq_t1e frequentist   P_delta=0(BF01 <= k1)
     freq_power_main frequentist P_delta=0.4(BF01 <= k1)
        estimate target direction       slack satisfied
     0.005608335   0.05        <= 0.044391665      TRUE
     0.803068443   0.80        >= 0.003068443      TRUE

We see that the number of patients decreases compared to the two-sided
test. We can plot the design:

``` r

plot(pulmonary_directional_freq)
```

![Figure 2: Output of the plot function for a one-sample Bayesian t-test
design based on frequentist
calibration.](figures/bfbin2arm-onesample-t-test-freq.png)

Figure 2: Output of the plot function for a one-sample Bayesian t-test
design based on frequentist calibration.

Figure 2 shows that the shape of the decision function now reflects the
one-sided test. For increasingly large t-statistic, the Bayes factor
will reject $`H_0`$. Note that the top right panel in Figure 2 shows
Bayesian power (as the design prior under $`H_1`$ is specified in the
common parameters above), but only frequentist power and type-I-error
are used for calibrating the required sample size.

### Bayesian calibration

This design controls the $`H_0`$-design-prior average probability of
incorrectly obtaining compelling evidence for benefit and requires 80%
average Bayesian power under the $`H_1`$ design prior.

``` r

pulmonary_directional_bayes <- do.call(
  calibrate_bf_ttest,
  c(
    directional_common_args,
    list(
      bayes_t1e = 0.05,
      bayes_power = 0.80
    )
  )
)

print(pulmonary_directional_bayes)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: greater.pointnull 
    n*: 134 (df = 133 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 2.666093   Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
      constraint     type
       bayes_t1e Bayesian
     bayes_power Bayesian
                    estimand   estimate target
     P_design,H0(BF01 <= k1) 0.00431281   0.05
     P_design,H1(BF01 <= k1) 0.80025788   0.80
     direction        slack satisfied
            <= 0.0456871901      TRUE
            >= 0.0002578833      TRUE

Since the Bayesian type-I-error under $`H_0:\delta = 0`$ is computed
under a point-prior, the Bayesian type-I-error rate is identical to the
frequentist one. This is the reason that the plot shows frequentist
type-I-error (which usually is more conservative), even though the
calibration uses Bayesian type-I-error (and power).

``` r

plot(pulmonary_directional_bayes)
```

![Figure 3: Output of the plot function for a one-sample Bayesian t-test
design based on Bayesian
calibration.](figures/bfbin2arm-onesample-t-test-bayes.png)

Figure 3: Output of the plot function for a one-sample Bayesian t-test
design based on Bayesian calibration.

### Hybrid calibration

This hybrid calibration uses a frequentist type-I error requirement at
$`\delta=0`$ and a Bayesian power requirement averaged over the $`H_1`$
design prior.

``` r

pulmonary_directional_bayes_power_freq_t1e <- do.call(
  calibrate_bf_ttest,
  c(
    directional_common_args,
    list(
      freq_t1e = 0.05,
      bayes_power = 0.80
    )
  )
)

print(pulmonary_directional_bayes_power_freq_t1e)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: greater.pointnull 
    n*: 134 (df = 133 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 2.666093   Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
      constraint        type
        freq_t1e frequentist
     bayes_power    Bayesian
                    estimand   estimate target
       P_delta=0(BF01 <= k1) 0.00431281   0.05
     P_design,H1(BF01 <= k1) 0.80025788   0.80
     direction        slack satisfied
            <= 0.0456871901      TRUE
            >= 0.0002578833      TRUE

### Full calibration without PCE

This version activates frequentist and Bayesian type-I error
requirements, as well as frequentist and Bayesian power requirements. It
does not impose a probability-of-compelling-evidence (PCE) requirement
for $`H_0`$.

``` r

pulmonary_directional_full_no_pce <- do.call(
  calibrate_bf_ttest,
  c(
    directional_common_args,
    list(
      freq_t1e = 0.05,
      bayes_t1e = 0.05,
      freq_power = list(
        main = list(effect = 0.40, target = 0.80)
      ),
      bayes_power = 0.80
    )
  )
)

print(pulmonary_directional_full_no_pce)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: greater.pointnull 
    n*: 134 (df = 133 )
    Thresholds: k1 = 0.1 , k0 = 3 
    Analysis prior under H0:
    Point prior: delta = 0
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 2.666093   Inf
    Evidence for H0 region(s) on t-scale:
    none

    Constraints:
          constraint        type
            freq_t1e frequentist
           bayes_t1e    Bayesian
         bayes_power    Bayesian
     freq_power_main frequentist
                    estimand   estimate target
       P_delta=0(BF01 <= k1) 0.00431281   0.05
     P_design,H0(BF01 <= k1) 0.00431281   0.05
     P_design,H1(BF01 <= k1) 0.80025788   0.80
     P_delta=0.4(BF01 <= k1) 0.97401497   0.80
     direction        slack satisfied
            <= 0.0456871901      TRUE
            <= 0.0456871901      TRUE
            >= 0.0002578833      TRUE
            >= 0.1740149710      TRUE

We see that all calibrations except the frequentit one require
$`n^{\star}=134`$ patients. Frequentist calibration under the much more
optimistic point prior at $`\delta = 0.4`$ requires solely
$`n^{\star}=75`$ patients.

The point-null versus greater directional test shows that the
directional Bayes factor is asymmetric in the observed $`t`$ statistic,
compare the top left panel in Figure 2. The region providing evidence
for $`H_1`$ should therefore occur in the positive tail only. The
implementation determines decision regions by finding the relevant
Bayes-factor roots and classifying the resulting intervals rather than
imposing a symmetric decision geometry.

## Example 3: Directional composite-null design

This example is inspired by an 8-week open-label, multicentre study of
inhaled tobramycin in people with cystic fibrosis and *Pseudomonas
aeruginosa* infection, compare Sutharsan et al. (2020). The published
study assessed pulmonary outcomes including forced expiratory volume in
1 second (FEV1) and lung-clearance index (LCI).

We consider a hypothetical future single-arm proof-of-concept study. Let
the endpoint be improvement in LCI from baseline to Week 8, coded as
baseline minus follow-up LCI. Thus, a positive value denotes clinical
improvement.

The scientific question is whether the intervention provides evidence of
benefit rather than merely evidence against an exact zero mean change:

``` math
H_0:\delta\leq0
\quad\text{versus}\quad
H_1:\delta>0.
```

The null-side hypothesis includes both no improvement and deterioration.
This is appropriate when a benefit claim should be made only when the
study yields compelling evidence for an improvement.

``` r

composite_lci_continuous <- calibrate_bf_ttest(
  design = "one.sample",
  hypothesis = "greater.composite",
  n_min = 20,
  n_max = 120,
  k1 = 1 / 20,
  k0 = 3,

  analysis_prior_h0 = make_analysis_prior_less_equal(sd = 0.5),
  analysis_prior_h1 = make_analysis_prior_greater(sd = 0.5),

  freq_t1e = 0.05,
  freq_power = list(main = list(effect = 0.40, target = 0.80)),
  bayes_t1e = 0.05,
  bayes_power = 0.80,
  pce_h0 = 0.60,

  design_prior_h0 = list(
    type = "truncnormal",
    mean = -0.10,
    sd = 0.25,
    lower = -Inf,
    upper = 0
  ),
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = fast_control
)

composite_lci_continuous
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: one.sample 
    Hypothesis: greater.composite 
    n*: 65 (df = 64 )
    Thresholds: k1 = 0.05 , k0 = 3 
    Analysis prior under H0:
    Prior: TN(0, 0.5^2; (-Inf, 0]) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, 0]
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 1.733227   Inf
    Evidence for H0 region(s) on t-scale:
         lower      upper
    [1,]  -Inf -0.6935961

    Constraints:
          constraint        type                estimand    estimate
            freq_t1e frequentist   P_delta=0(BF01 <= k1) 0.043935585
           bayes_t1e    Bayesian P_design,H0(BF01 <= k1) 0.005294749
         bayes_power    Bayesian P_design,H1(BF01 <= k1) 0.802264572
              pce_h0    Bayesian P_design,H0(BF01 >= k0) 0.758721651
     freq_power_main frequentist P_delta=0.4(BF01 <= k1) 0.930717800
     target direction       slack satisfied
       0.05        <= 0.006064415      TRUE
       0.05        <= 0.044705251      TRUE
       0.80        >= 0.002264572      TRUE
       0.60        >= 0.158721651      TRUE
       0.80        >= 0.130717800      TRUE

``` r

plot(composite_lci_continuous)
```

![Figure 3: Output of the plot function for a one-sample Bayesian
composite-null versus directional t-test design based on full
calibration. Both Bayesian and frequentist power and type-I-error are
calibrated jointly. Also, a 60% threshold on the probability of
compelling evidence for the null hypothesis is
guaranteed.](figures/bfbin2am-onesample-t-test-ex3.png)

Figure 3: Output of the plot function for a one-sample Bayesian
composite-null versus directional t-test design based on full
calibration. Both Bayesian and frequentist power and type-I-error are
calibrated jointly. Also, a 60% threshold on the probability of
compelling evidence for the null hypothesis is guaranteed.

As a sidenote, single-arm change-from-baseline evidence is vulnerable to
within-person variability, natural history, regression to the mean, and
measurement error; it does not itself establish a causal treatment
effect, in contrast to two-arm change-from-baseline evidence. However,
the composite-null formulation clarifies the evidential decision rule.
What it cannot achieve, however, it to remove the limitations of a
non-comparative design. External controls, randomization, or other
design features may be necessary for a confirmatory efficacy claim. We
study this example nevertheless, because a directional composite-null
design is widely used in clinical trials.

## Comparison with a conventional frequentist design

The conventional comparator is a one-sided one-sample t-test,

``` math
H_0: \delta \leq 0
\qquad \text{versus} \qquad
H_1: \delta > 0,
```

designed for a type-I error rate of 0.05 and 80% power at the fixed
alternative $`\delta = 0.40`$. We also consider a Bayes-factor design
that retains the frequentist type-I-error and power requirements but
does not require compelling evidence for the H0 region. This separates
the sample-size implication of requiring strong evidence for H1 from the
additional requirement to obtain compelling H0 evidence.

``` r

## Conventional frequentist one-sided one-sample t-test
freq_design <- power.t.test(
  delta = 0.40,
  sd = 1,
  sig.level = 0.05,
  power = 0.80,
  type = "one.sample",
  alternative = "one.sided"
)

freq_n <- ceiling(freq_design$n)
```

``` r

## Bayes-factor design: H1 evidence only, no PCE requirement.
composite_lci_no_pce <- calibrate_bf_ttest(
  design = "one.sample",
  hypothesis = "greater.composite",
  n_min = 20,
  n_max = 120,
  k1 = 1 / 20,
  k0 = 3,

  analysis_prior_h0 = make_analysis_prior_less_equal(sd = 0.5),
  analysis_prior_h1 = make_analysis_prior_greater(sd = 0.5),

  freq_t1e = 0.05,
  freq_power = list(main = list(effect = 0.40, target = 0.80)),

  design_prior_h0 = list(
    type = "truncnormal",
    mean = -0.10,
    sd = 0.25,
    lower = -Inf,
    upper = 0
  ),
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = fast_control
)
```

``` r

sample_size_comparison <- data.frame(
  design = c(
    "Frequentist one-sample t-test",
    "Bayes-factor design: H1 evidence only",
    "Bayes-factor design: H1 and H0 evidence"
  ),
  type_I_error_control = c(
    "One-sided alpha = 0.05",
    "P(delta = 0; BF01 <= 1/20) <= 0.05",
    "P(delta = 0; BF01 <= 1/20) <= 0.05"
  ),
  power_requirement = c(
    "Power(delta = 0.40) >= 0.80",
    "P(delta = 0.40; BF01 <= 1/20) >= 0.80",
    "P(delta = 0.40; BF01 <= 1/20) >= 0.80"
  ),
  H0_evidence_requirement = c(
    "None",
    "None",
    "P_design,H0(BF01 >= 3) >= 0.60"
  ),
  sample_size = c(
    freq_n,
    composite_lci_no_pce$n,
    composite_lci_continuous$n
  ),
  stringsAsFactors = FALSE
)

knitr::kable(
  sample_size_comparison,
  col.names = c(
    "Design",
    "Type-I-error control",
    "Power requirement",
    "H0-evidence requirement",
    "Required n"
  ),
  caption = "Sample-size comparison under a one-sided standardized effect of 0.40."
)
```

| Design | Type-I-error control | Power requirement | H0-evidence requirement | Required n |
|:---|:---|:---|:---|---:|
| Frequentist one-sample t-test | One-sided alpha = 0.05 | Power(delta = 0.40) \>= 0.80 | None | 41 |
| Bayes-factor design: H1 evidence only | P(delta = 0; BF01 \<= 1/20) \<= 0.05 | P(delta = 0.40; BF01 \<= 1/20) \>= 0.80 | None | 43 |
| Bayes-factor design: H1 and H0 evidence | P(delta = 0; BF01 \<= 1/20) \<= 0.05 | P(delta = 0.40; BF01 \<= 1/20) \>= 0.80 | P_design,H0(BF01 \>= 3) \>= 0.60 | 65 |

Sample-size comparison under a one-sided standardized effect of 0.40.
{.table style="width:100%;"}

The conventional frequentist calculation requires 41 patients. The BF
design without a probability-of-compelling-evidence (PCE) requirement
requires 43 patients. Its slightly larger sample size reflects the use
of the more stringent evidential criterion $`BF_{01} \leq 1/20`$, rather
than a conventional one-sided significance threshold.

Adding the requirement that the design yield compelling $`H_0`$ evidence
($`BF_{01} \geq 3`$) with probability at least 0.60 under the $`H_0`$
design prior increases the required sample size to 65. Thus, the PCE
requirement adds 22 patients beyond the $`H_1`$-evidence-only BF design,
and the complete calibrated design adds 24 patients beyond the
conventional frequentist design.

## Example 4: Paired one-sample design

This example is inspired by the Healthy BP study, which compared a
dietitian-led telehealth lifestyle intervention with a self-directed
intervention. Its primary endpoint was change in 24-hour ambulatory
systolic blood pressure from baseline to 12 weeks, see Taher et al.
(2021).

Here, we consider a hypothetical future single-group implementation
study with measurements at the same two time points. This is a
prospective design illustration, not a reproduction of the original
randomized trial or its between-group analysis. All effect sizes, prior
distributions, evidence thresholds, and calibration targets below are
illustrative planning assumptions rather than quantities attributed to
that trial.

### Paired outcome and hypothesis

Define the participant-level outcome as

``` math
  D_i
  =
  X_{i,\mathrm{baseline}}
  -
  X_{i,\mathrm{12\ weeks}},
```

so that a positive difference represents a beneficial reduction in
systolic blood pressure. Let

``` math
  \mu_D = \mathbb{E}(D_i),
  \qquad
  \sigma_D = \operatorname{SD}(D_i),
  \qquad
  \delta = \frac{\mu_D}{\sigma_D}.
```

A clinically relevant reduction of 5 mmHg and an assumed planning
standard deviation of paired differences of 12.5 mmHg imply

``` math
  \delta_A = \frac{5}{12.5} = 0.40.
```

The standard deviation in this calculation is the standard deviation of
the within-participant differences, not the standard deviation of
baseline or follow-up blood pressure considered separately.

We use the directional composite-null hypotheses

``` math
  H_0:\delta\leq0
  \qquad\text{versus}\qquad
  H_1:\delta>0.
```

Thus, the null-side hypothesis includes both no mean reduction and a
mean increase in blood pressure.

For a paired design, $`n`$ is the number of participants with evaluable
baseline and follow-up measurements. The test statistic is

``` math
  T = \frac{\sqrt{n}\,\bar D}{S_D},
  \qquad
  T\mid\delta \sim t_{n-1}(\sqrt{n}\delta).
```

The statistical engine is therefore the same as for a one-sample test of
the participant-level differences.

### Priors, evidence thresholds, and calibration targets

We use normal analysis priors with scale 0.5, truncated and normalized
on the respective hypothesis regions:

``` math
  \delta\mid H_0
  \sim N(0,0.5^2)\big|_{(-\infty,0]},
  \qquad
  \delta\mid H_1
  \sim N(0,0.5^2)\big|_{(0,\infty)}.
```

These analysis priors define the Bayes factor. For planning, we use
separate design priors:

``` math
  \pi_{d,0}
  =
  N(-0.10,0.25^2)\big|_{(-\infty,0]},
  \qquad
  \pi_{d,1}
  =
  N(0.40,0.20^2)\big|_{(0,\infty)}.
```

The evidence thresholds are fixed in advance:

- $`BF_{01}\leq1/20`$: compelling evidence for a beneficial reduction;
- $`BF_{01}\geq3`$: compelling evidence for the null-side hypothesis;
- otherwise: an inconclusive result.

We jointly require:

- pointwise frequentist type-I error at $`\delta=0`$ at most 0.05;
- frequentist power at $`\delta=0.40`$ at least 0.80;
- Bayesian type-I error under the null-side design prior at most 0.05;
- Bayesian power under the beneficial-effect design prior at least 0.80;
- probability of compelling null-side evidence under the null-side
  design prior at least 0.60.

``` r

paired_bp <- calibrate_bf_ttest(
  design = "paired",
  hypothesis = "greater.composite",

  # Number of evaluable participant pairs.
  n_min = 20,
  n_max = 120,

  # Prespecified, fixed evidence thresholds.
  k1 = 1 / 20,
  k0 = 3,

  # Default composite-null analysis priors:
  # N(0, 0.5^2), truncated to the respective half-lines.
  analysis_prior_sd = 0.5,

  freq_t1e = 0.05,
  freq_power = list(
    main = list(effect = 0.40, target = 0.80)
  ),
  bayes_t1e = 0.05,
  bayes_power = 0.80,
  pce_h0 = 0.60,

  design_prior_h0 = list(
    type = "truncnormal",
    mean = -0.10,
    sd = 0.25,
    lower = -Inf,
    upper = 0
  ),
  design_prior_h1 = list(
    type = "truncnormal",
    mean = 0.40,
    sd = 0.20,
    lower = 0,
    upper = Inf
  ),

  control = control_full_fast
)

print(paired_bp)
```

    Calibrated Bayesian t-test design
    ---------------------------------
    Design: paired 
    Hypothesis: greater.composite 
    n*: 65 (df = 64 )
    Thresholds: k1 = 0.05 , k0 = 3 
    Analysis prior under H0:
    Prior: TN(0, 0.5^2; (-Inf, 0]) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [-Inf, 0]
    Analysis prior under H1:
    Prior: TN(0, 0.5^2; [0, Inf)) 
      family:   normal
      location: 0
      scale:    0.5
      support:  [0, Inf]
    Evidence for H1 region(s) on t-scale:
            lower upper
    [1,] 1.733227   Inf
    Evidence for H0 region(s) on t-scale:
         lower      upper
    [1,]  -Inf -0.6935961

    Constraints:
          constraint        type                estimand    estimate
            freq_t1e frequentist   P_delta=0(BF01 <= k1) 0.043935585
           bayes_t1e    Bayesian P_design,H0(BF01 <= k1) 0.005294749
         bayes_power    Bayesian P_design,H1(BF01 <= k1) 0.802264572
              pce_h0    Bayesian P_design,H0(BF01 >= k0) 0.758721651
     freq_power_main frequentist P_delta=0.4(BF01 <= k1) 0.930717800
     target direction       slack satisfied
       0.05        <= 0.006064415      TRUE
       0.05        <= 0.044705251      TRUE
       0.80        >= 0.002264572      TRUE
       0.60        >= 0.158721651      TRUE
       0.80        >= 0.130717800      TRUE

### Calibrated design

The selected sample size is $`n^\star = 65`$ evaluable participant
pairs, with $`64`$ degrees of freedom. Thus, the design requires at
least $`65`$ participants contributing both measurements, not
$`2\cdot 65=130`$ independently sampled participants.

### Interpretation and relation to the one-sample example

This example uses the same standardized-effect priors, fixed evidence
thresholds, and calibration requirements as the composite-null
one-sample example. Consequently, with identical numerical controls, the
paired and one-sample calculations represent the same mathematical
calibration problem. What changes is the interpretation of the outcome
and sample size: here the observations are participant-level differences
and $`n`$ counts evaluable pairs.

Repeated baseline and follow-up measurements are not independent
observations. Their dependence is incorporated through the standard
deviation of the paired differences. If planning information is
available only for the marginal standard deviations and
within-participant correlation, then

``` math
  \sigma_D^2
  =
  \sigma_{\mathrm{baseline}}^2
  +
  \sigma_{\mathrm{follow\text{-}up}}^2
  -
  2\rho\,
  \sigma_{\mathrm{baseline}}
  \sigma_{\mathrm{follow\text{-}up}}.
```

The assumed value $`\sigma_D=12.5`$ mmHg should therefore be justified
using relevant paired data or examined in sensitivity analyses.

Finally, the calibrated sample size refers to evaluable pairs.
Recruitment inflation for missing follow-up measurements should be
considered separately. Neither such inflation nor the paired analysis
itself resolves bias arising from informative missingness or the causal
limitations of a single-group study.

## Interpretation

[`calibrate_bf_ttest()`](https://rikokelter.github.io/bfbin2arm/reference/calibrate_bf_ttest.md)
returns the smallest sample size in the requested range that satisfies
the activated constraints. The returned object includes the Bayes-factor
evidence regions, attained operating characteristics, constraint table,
and numerical controls used for the calculation.

The examples are illustrations of a design procedure. They do not imply
that a Bayesian design automatically requires fewer participants than a
conventional design; any difference follows from the chosen evidence
thresholds and jointly calibrated operating-characteristic requirements.

## References

Sutharsan, Sivagurunathan, Susanne Naehrig, Uwe Mellies, Christian
Sieder, and Jörg Ziegler. 2020. “An 8 Week Open-Label Interventional
Multicenter Study to Explore the Lung Clearance Index as Endpoint for
Clinical Trials in Cystic Fibrosis Patients $`\geq`$8 Years of Age,
Chronically Infected with Pseudomonas Aeruginosa.” *BMC Pulm Med* 20:
167. <https://doi.org/10.1186/s12890-020-01201-y>.

Taher, Mohamed, Christina Yule, Heather Bonaparte, et al. 2021.
“Telehealth Versus Self-Directed Lifestyle Intervention to Promote
Healthy Blood Pressure: A Protocol for a Randomised Controlled Trial.”
*BMJ Open* 11 (3): e044292.
<https://doi.org/10.1136/bmjopen-2020-044292>.

Vijaykumar, Kadambari, George M. Solomon, Jennifer Guimbellot, et al.
2025. “Ivacaftor for Chronic Obstructive Pulmonary Disease: Results from
a Phase 2, Randomized Controlled Trial.” *Am J Respir Crit Care Med* 211
(5): 823–31. <https://doi.org/10.1164/rccm.202407-1302OC>.
