# ROPE-based trial design for single-arm one-stage phase II trials with binary endpoints

## Introduction

This vignette illustrates how to use
[`design_singlearm_onestage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_onestage_rope.md)
to calibrate ROPE-based equivalence designs for single-arm phase II
trials with binary endpoints. ROPE stands for the region of practical
equivalence and has been proposed by Kruschke and Liddell (2018),
Kruschke (2014) and Kruschke (2018), even though the idea itself is
older and appears under various names in different contexts, see Kelter
(2021), Liao et al. (2020), Linde et al. (2023), Lakens et al. (2018),
Wellek (2010) and Pan et al. (2025). The idea to replace the test of a
point-null hypothesis with a small-interval goes at least back until
Hodges and Lehmann (1954).

## Setup

We consider a single-arm binomial model

``` math
Y \mid p \sim \mathrm{Binomial}(n, p),
```

where $`Y`$ is the number of responders among $`n`$ patients and
$`p \in (0,1)`$ is the true response probability under the experimental
treatment. We fix a benchmark response rate $`p_0`$ (e.g. historical
control or standard of care) and define the risk difference

``` math
\Delta = p - p_0
```

We work with a symmetric ROPE formulation on the risk-difference scale.
Let $`\Delta = p - p_0`$ denote the risk difference between the
experimental treatment and the benchmark response probability $`p_0`$,
and let $`\delta > 0`$ be the equivalence margin. On the risk-difference
scale we define

``` math
H_0:\; |\Delta| > \delta,
```

``` math
H_1:\; |\Delta| \le \delta.
```

Equivalently, on the response-probability scale the ROPE is

``` math
[p_0 - \delta,\; p_0 + \delta] \cap (0,1),
```

and the hypotheses can be written as

``` math
H_0:\; p \notin [\,p_0 - \delta,\; p_0 + \delta\,],
```

``` math
H_1:\; p \in [\,p_0 - \delta,\; p_0 + \delta\,].
```

## The region of practical equivalence (ROPE)

The **region of practical equivalence (ROPE)** on the risk-difference
scale is

``` math
\mathcal{R}_\Delta = [-\delta, \delta],
```

where $`\delta > 0`$ is a prespecified equivalence margin. Equivalently,
on the response-probability scale the ROPE for $`p`$ is

``` math
\mathcal{R}_p = [p_0 - \delta,\; p_0 + \delta] \cap (0,1).
```

Given a beta analysis prior

``` math
p \sim \mathrm{Beta}(a, b),
```

the posterior after observing $`Y = y`$ is

``` math
p \mid y \sim \mathrm{Beta}(a + y,\; b + n - y),
```

and the posterior ROPE probability is

``` math
\Pr\bigl(p \in \mathcal{R}_p \mid y\bigr)
  = F_{\mathrm{Beta}(a+y,\,b+n-y)}(p_0 + \delta)
  - F_{\mathrm{Beta}(a+y,\,b+n-y)}(p_0 - \delta),
```

with endpoints truncated to $`[0,1]`$ if needed. A ROPE-based
equivalence decision rule declares **practical equivalence** if

``` math
\Pr\bigl(p \in \mathcal{R}_p \mid y\bigr) \ge \gamma_{\mathrm{eq}},
```

where $`\gamma_{\mathrm{eq}} \in (0.5, 1)`$ is a pre-specified evidence
threshold.

## Design and analysis priors

At the **design stage** we distinguish between three priors:

- an **analysis prior** $`\mathrm{Beta}(a, b)`$ used to compute
  posterior ROPE probabilities,
- a **design prior under equivalence**
  $`H_1: \Delta \in [-\delta, \delta]`$, typically
  $`\mathrm{Beta}(a_1, b_1)`$ centred near $`p_0`$,
- a **design prior under non-equivalence**
  $`H_0: \Delta \notin [-\delta, \delta]`$, typically
  $`\mathrm{Beta}(a_0, b_0)`$ centred away from the ROPE.

These design priors induce beta–binomial predictive distributions for
$`Y`$ under equivalence and non-equivalence, respectively. Under the
equivalence design prior $`\pi_1`$ we define **ROPE-based Bayesian
power** as

``` math
\text{Power}_\text{ROPE}(n)
= \Pr_{\pi_1}\bigl( \Pr(p \in \mathcal{R}_p \mid Y) \ge \gamma_{\mathrm{eq}} \bigr),
```

and under the non-equivalence design prior $`\pi_0`$ we define the
**ROPE-based Bayesian type-I error** as

``` math
\alpha_\text{ROPE}(n)
= \Pr_{\pi_0}\bigl( \Pr(p \in \mathcal{R}_p \mid Y) \ge \gamma_{\mathrm{eq}} \bigr).
```

## ROPE decision illustrations

In this section we illustrate the ROPE-based decision rule for four
prototypical outcomes in a single-arm binomial model with analysis prior
$`p \sim \mathrm{Beta}(1,1)`$, benchmark response rate $`p_0 = 0.30`$,
and ROPE $`\mathcal{R}_p = [p_0 - \delta, p_0 + \delta] = [0.18, 0.42]`$
with $`\delta = 0.12`$.

For an observed responder count $`Y = y`$ out of $`n`$ patients, the
posterior is

``` math
p \mid y \sim \mathrm{Beta}(a + y,\; b + n - y),
```

and the symmetric ROPE probability is

``` math
\Pr\bigl(|p - p_0| \le \delta \mid y\bigr)
  = \Pr(p_0 - \delta \le p \le p_0 + \delta \mid y).
```

We adopt the following simple decision rule:

- **Equivalence accepted** if
  $`\Pr(|p - p_0| \le \delta \mid y) \ge \gamma_{\mathrm{eq}}`$.
- **Non-equivalence accepted** if
  $`\Pr(|p - p_0| > \delta \mid y) \ge \gamma_{\mathrm{diff}}`$.
- **Indecisive** otherwise,

with $`\gamma_{\mathrm{eq}} = \gamma_{\mathrm{diff}} = 0.80`$ in the
examples below.

### 1) Equivalence accepted

We choose an outcome $`(n, y)`$ for which the posterior is concentrated
inside the ROPE and
$`\Pr(|p - p_0| \le \delta \mid y) \ge \gamma_{\mathrm{eq}}`$, so the
decision is to **accept equivalence**.

![Figure 1: Illustration of the first possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Equivalence is accepted,
because sufficient posterior probability mass concentrates inside the
ROPE. The true data-generating process follows the alternative
hypothesis, that is, equivalence indeed
holds.](figures/singlearm-onestage-rope-scenario1.png)

Figure 1: Illustration of the first possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Equivalence is accepted,
because sufficient posterior probability mass concentrates inside the
ROPE. The true data-generating process follows the alternative
hypothesis, that is, equivalence indeed holds.

The plot illustrates this first possible outcome.

### 2) Type-I error: equivalence concluded under the null hypothesis

Conceptually, a type-I error occurs when the *true* data-generating
process is non-equivalent (e.g. $`p = 0.55`$ or 0.60), but the observed
data still lead the ROPE rule to **accept equivalence**. Thus, $`H_0`$
is true and $`p \notin [\,p_0 - \delta,\; p_0 + \delta\,]`$ holds.

In this plot we **do not** change the posterior calculation—posterior is
always conditional on the observed $`(n,y)`$ and the analysis prior. To
illustrate a type-I error, we choose $`(n,y)`$ such that:

- $`y`$ is plausible under a non-equivalence scenario (e.g. generated
  from $`p = 0.55`$), **and**
- the resulting posterior still satisfies
  $`\Pr(|p - p_0| \le \delta \mid y) \ge \gamma_{\mathrm{eq}}`$.

For illustration we tune $`y`$ so that this happens:

![Figure 2: Illustration of the second possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Equivalence is accepted,
because sufficient posterior probability mass concentrates inside the
ROPE. In contrast to the first possible scenario, the true
data-generating process follows the null hypothesis. Thus, a ROPE-based
type-I-error occurs.](figures/singlearm-onestage-rope-scenario2.png)

Figure 2: Illustration of the second possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Equivalence is accepted,
because sufficient posterior probability mass concentrates inside the
ROPE. In contrast to the first possible scenario, the true
data-generating process follows the null hypothesis. Thus, a ROPE-based
type-I-error occurs.

In this scenario, in contrast to scenario 1 above, the *true* $`p`$ lies
outside the ROPE (under $`H_0`$), but due to sampling variability the
posterior still concentrates enough mass inside the ROPE to meet the
equivalence threshold.

### 3) Indecisive result

Here we choose $`(n,y)`$ such that neither threshold is reached:

- $`\Pr(|p - p_0| \le \delta \mid y) < \gamma_{\mathrm{eq}}`$,
- $`\Pr(|p - p_0| > \delta \mid y) < \gamma_{\mathrm{diff}}`$.

The posterior spreads substantial mass both inside and outside the ROPE,
and the decision is **indecisive**.

![Figure 3: Illustration of the third possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: The result is indecisive,
because neither does sufficient posterior probability mass concentrate
inside the ROPE, nor outside the
ROPE.](figures/singlearm-onestage-rope-scenario3.png)

Figure 3: Illustration of the third possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: The result is indecisive,
because neither does sufficient posterior probability mass concentrate
inside the ROPE, nor outside the ROPE.

### 4) Clear non-equivalence

Finally, we choose an outcome where the posterior lies mostly outside
the ROPE, so that
$`\Pr(|p - p_0| > \delta \mid y) \ge \gamma_{\mathrm{diff}}`$ and we
**accept non-equivalence**.

![Figure 4: Illustration of the fourth possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Non-equivalence is
accepted, because sufficient posterior probability mass concentrates
outside the ROPE.](figures/singlearm-onestage-rope-scenario4.png)

Figure 4: Illustration of the fourth possible scenario in a ROPE-based
clinical phase II trial with binary endpoints: Non-equivalence is
accepted, because sufficient posterior probability mass concentrates
outside the ROPE.

In this last case, the treatment is worse than the standard of care with
success probability $`p_0`$.

## A first ROPE-based design example

We now provide a simple example of the calibration function
[`design_singlearm_onestage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_onestage_rope.md),
which calibrates a single-arm one-stage phase II design using the ROPE
as the primary measure of evidence.

We consider a setting with benchmark response rate $`p_0 = 0.30`$ and
regard differences up to 0.12 as clinically negligible. Thus the ROPE on
$`p`$ is $`\mathcal{R}_p = [0.18, 0.42]`$.

We use:

- a **uniform analysis prior** $`\mathrm{Beta}(1,1)`$,
- a **non-equivalence design prior** $`\mathrm{Beta}(60,40)`$ with mean
  0.60, representing clearly superior response compared to 0.30
  (non-equivalence); this is the design prior under $`H_0`$
- an **equivalence design prior** $`\mathrm{Beta}(36,84)`$ with mean
  0.30, representing plausible equivalence scenarios; this is the design
  prior under $`H_1`$
- an **equivalence threshold** $`\gamma_{\mathrm{eq}} = 0.80`$,
- a **target ROPE-based power** of 0.80 under the equivalence design
  prior,
- a **maximum ROPE-based type-I error** of 0.10 under the
  non-equivalence design prior,
- a **sustain requirement** of `sustain_n = 10`, meaning the criteria
  must hold for 10 consecutive sample sizes starting from the selected
  $`n^\ast`$.

``` r

des_baseline <- design_singlearm_onestage_rope(
  n_min = 20,
  n_max = 200,
  p0 = 0.30,     # benchmark response rate p0
  delta = 0.12,  # ROPE half-width: equivalence if p in [0.18, 0.42]
  gamma_eq = 0.80,  # posterior ROPE probability threshold for equivalence

  # Analysis prior: p ~ Beta(a, b), used for posterior and ROPE decision
  a = 1,
  b = 1,

  # Design prior under H0 (non-equivalence): p ~ Beta(da0, db0)
  # Here: mean 0.60, representing clearly higher response than 0.30.
  da0 = 60,
  db0 = 40,

  # Design prior under H1 (equivalence): p ~ Beta(da1, db1)
  # Here: mean 0.30, representing plausible equivalence scenarios.
  da1 = 36,
  db1 = 84,

  # Target ROPE-based power under H1 (equivalence design prior)
  target_power = 0.80,

  # Maximum ROPE-based type-I error under H0 (non-equivalence design prior)
  target_type1 = 0.10,

  # Stability requirement: criteria must hold for 10 consecutive n values
  sustain_n = 10
)
```

We can take a look at the resulting design object:

``` r

des_baseline
#> One-stage single-arm ROPE design
#> Direction: equivalence 
#> Calibration: Bayesian 
#> Search range n: 20 to 200 
#> Null probability p0: 0.3 
#> Margin delta: 0.12 
#> Probability threshold gamma_eq: 0.8 
#> Probability threshold gamma_diff: 0.8 
#> Analysis prior: Beta(1, 1)
#> Design prior (H0): Beta(60, 40)
#> Design prior (H1): Beta(36, 84)
#> Target Bayesian power: 0.8 
#> Target Bayesian type-I error: 0.1 
#> Sustain n: 10 
#> Selected sample size n*: 94 
#> Bayesian power(n*): 0.8231 
#> Bayesian type-I(n*): 0.0009 
#> PCE(H0)(n*): 0.9730 
#> Equivalence region: {20-35}
#> Compelling evidence for non-equivalence region: {0-13, 44-94}
```

The printed output reports:

- the search range for $`n`$,
- the ROPE specification (`p0`, `delta`, `gamma_eq`),
- the analysis and design priors in beta parameterization,
- the target power and type-I constraints,
- the chosen `sustain_n`,
- the selected sample size `Selected n`,
- the ROPE-based power and type-I error at that $`n`$,
- and the equivalence decision region
  $`[y_{\min}^{\mathrm{eq}}, y_{\max}^{\mathrm{eq}}]`$, i.e. all
  responder counts $`y`$ that lead to practical equivalence.

### Summarizing the design

We can summarize the calibration grid and the selected design via:

``` r

summary(des_baseline)
```

The summary object (not shown here) contains:

- the selected row of the grid (with `n`, `y_eq_min`, `y_eq_max`,
  `power`, `type1`),
- the first and last 10 rows of the evaluated `n` values.

In particular:

- `y_eq_min` and `y_eq_max` are the smallest and largest responder
  counts for which the posterior ROPE probability exceeds `gamma_eq` and
  equivalence would be concluded;
- `power` is the ROPE-based Bayesian power under the equivalence design
  prior at that `n`;
- `type1` is the ROPE-based Bayesian type-I error under the
  non-equivalence design prior at that `n`.

These summaries allow you to inspect how power and type-I error evolve
with increasing sample size, and how the equivalence decision region
moves.

### Plotting the design

The overview plot visualizes operating characteristics, priors, and a
textual summary:

``` r

plot(des_baseline)
```

![Figure 5: Illustration of calibrated single-arm one-stage design of a
ROPE-based clinical phase II trial with binary
endpoint.](figures/singlearm-onestage-rope-fig5.png)

Figure 5: Illustration of calibrated single-arm one-stage design of a
ROPE-based clinical phase II trial with binary endpoint.

- The **upper left panel** shows ROPE-based power and type-I error as
  functions of $`n`$, with horizontal lines at `target_power` and
  `target_type1`, and a vertical line at the selected `n`.
- The **upper right panel** displays a textual summary of the key inputs
  and outputs (priors, ROPE, thresholds, selected `n`, power, type-I,
  and equivalence region).
- The **lower left panel** displays the design priors under `H0` and
  `H1` overlaid: their beta densities highlight which response
  probabilities are regarded as typical under non-equivalence and
  equivalence, respectively.
- The **lower right panel** displays the analysis prior
  $`\mathrm{Beta}(a,b)`$, which governs the posterior ROPE probabilities
  used in the decision rule.

You can also visualize only the operating characteristics or the
decision region:

``` r

plot(des_baseline, what = "operating_characteristics")
```

![Figure 6: Visualization of the operating characteristics of a
calibrated single-arm one-stage design of a ROPE-based clinical phase II
trial with binary endpoint.](figures/singlearm-onestage-rope-fig6.png)

Figure 6: Visualization of the operating characteristics of a calibrated
single-arm one-stage design of a ROPE-based clinical phase II trial with
binary endpoint.

``` r

plot(des_baseline, what = "decision_region")
```

![Figure 7: Visualization of the equivalence region for increasing
sample size of a calibrated single-arm one-stage design of a ROPE-based
clinical phase II trial with binary
endpoint.](figures/singlearm-onestage-rope-fig7.png)

Figure 7: Visualization of the equivalence region for increasing sample
size of a calibrated single-arm one-stage design of a ROPE-based
clinical phase II trial with binary endpoint.

The decision-region plot shows how the range of responder counts leading
to equivalence changes with `n`, providing intuition about how stringent
the rule is at different sample sizes.

## Example 1: Oncology phase II equivalence trial

In this section we illustrate a full ROPE-based design calibration in a
setting resembling a single-arm phase II oncology trial with a binary
endpoint such as objective response rate (ORR), compare Chen et al.
(2022), Kelter and Schnurr (2024) and Lee and Liu (2008). For
definiteness, we assume:

- Historical control ORR $`p_0 = 0.25`$ based on previous phase II data.
- The new treatment is considered *clinically non-inferior / equivalent*
  if its true ORR lies within ±12 percentage points of $`p_0`$, that is,
  $`\mathcal{R}_p = [0.13, 0.37]`$. This is a common margin in phase II
  oncology trials, compare Hashim et al. (2021).
- We want a high probability to conclude practical equivalence when the
  true ORR is near 0.25, and a low probability to conclude equivalence
  when the true ORR is clearly better or worse than 0.25
  (non-equivalence).

### Clinical hypotheses and ROPE

On the response-probability scale we set $`p_0 = 0.30`$ and
$`\delta = 0.12`$. The ROPE for equivalence is

``` math
\mathcal{R}_p = [p_0 - \delta,\; p_0 + \delta] = [0.18, 0.42].
```

We formulate the hypotheses as

``` math
H_0:\; p \notin [\,p_0 - \delta,\; p_0 + \delta\,]
     \quad \text{(non-equivalence, clinically relevant difference)},
```
``` math
H_1:\; p \in [\,p_0 - \delta,\; p_0 + \delta\,]
     \quad \text{(practical equivalence)}.
```

We adopt the following ROPE-based decision rule:

- **Accept equivalence** ($`H_1`$) if
  $`\Pr(p \in \mathcal{R}_p \mid y) \ge \gamma_{\mathrm{eq}}`$.
- **Accept non-equivalence** ($`H_0`$) if
  $`\Pr(p \notin \mathcal{R}_p \mid y) \ge \gamma_{\mathrm{diff}}`$.
- **Indecisive** otherwise.

For this example, we set
$`\gamma_{\mathrm{eq}} = \gamma_{\mathrm{diff}} = 0.80`$.

### Analysis and design priors

We separate the analysis prior from the design priors.

- **Analysis prior** for ORR:

  ``` math
  p \sim \mathrm{Beta}(1,1),
  ```

  a uniform prior on $`(0,1)`$, reflecting weak prior information.

- **Design prior under equivalence** $`H_1`$:

  ``` math
  p \sim \mathrm{Beta}(a_1, b_1) = \mathrm{Beta}(36, 84),
  ```

  which has mean $`36 / (36 + 84) = 0.30`$ and moderate concentration
  around $`p_0 = 0.30`$. This prior represents plausible ORR values
  under practical equivalence.

- **Design prior under non-equivalence** $`H_0`$: we consider superior
  scenarios where ORR is clinically higher than 0.42. For concreteness
  we choose

  ``` math
  p \sim \mathrm{Beta}(60, 40),
  ```

  which is centred at 0.6 and places most mass clearly outside the ROPE
  interval \[0.18, 0.42\]. This prior represents clinically relevant
  departures from equivalence (e.g. strong improvement), and is used to
  quantify ROPE-based type-I error for wrongly declaring equivalence in
  such scenarios.

These design priors induce beta–binomial predictive distributions for
the response count $`Y`$ under $`H_1`$ and $`H_0`$, respectively.

Under the equivalence design prior $`\pi_1`$, the ROPE-based Bayesian
power is

``` math
\text{Power}_\text{ROPE}(n)
= \Pr_{\pi_1}\bigl( \Pr(p \in \mathcal{R}_p \mid Y) \ge \gamma_{\mathrm{eq}} \bigr),
```

and under the non-equivalence design prior $`\pi_0`$, the ROPE-based
Bayesian type-I error is

``` math
\alpha_\text{ROPE}(n)
= \Pr_{\pi_0}\bigl( \Pr(p \in \mathcal{R}_p \mid Y) \ge \gamma_{\mathrm{eq}} \bigr).
```

### Calibration target

For this oncology-inspired example we consider the following calibration
goals:

- ROPE-based power under $`H_1`$ at least 80%:
  $`\text{Power}_\text{ROPE}(n) \ge 0.80`$.
- ROPE-based type-I error under $`H_0`$ at most 10%:
  $`\alpha_\text{ROPE}(n) \le 0.10`$.
- A stability requirement `sustain_n = 10`, meaning that the criteria
  must hold for 10 consecutive sample sizes starting at the selected
  $`n^\ast`$. This guards against local non-monotonicities in the
  discrete predictive curves.

We search over a one-stage sample size range of 20 to 200 patients.

``` r

des_onc <- design_singlearm_onestage_rope(
  n_min = 20,
  n_max = 200,
  p0 = 0.30,
  delta = 0.12,
  gamma_eq = 0.80,

  # Analysis prior p ~ Beta(a, b)
  a = 1, b = 1,

  # Design priors under H0 and H1
  da0 = 60, db0 = 40,   # H0: non-equivalence, mean ~0.60
  da1 = 36, db1 = 84,   # H1: equivalence, mean ~0.3

  target_power = 0.80,
  target_type1 = 0.10,
  sustain_n = 10
)

des_onc
#> One-stage single-arm ROPE design
#> Direction: equivalence 
#> Calibration: Bayesian 
#> Search range n: 20 to 200 
#> Null probability p0: 0.3 
#> Margin delta: 0.12 
#> Probability threshold gamma_eq: 0.8 
#> Probability threshold gamma_diff: 0.8 
#> Analysis prior: Beta(1, 1)
#> Design prior (H0): Beta(60, 40)
#> Design prior (H1): Beta(36, 84)
#> Target Bayesian power: 0.8 
#> Target Bayesian type-I error: 0.1 
#> Sustain n: 10 
#> Selected sample size n*: 94 
#> Bayesian power(n*): 0.8231 
#> Bayesian type-I(n*): 0.0009 
#> PCE(H0)(n*): 0.9730 
#> Equivalence region: {20-35}
#> Compelling evidence for non-equivalence region: {0-13, 44-94}
```

The printed output shows the selected sample size $`n^\ast`$, ROPE-based
power and type-I error at that $`n^\ast`$, and the equivalence decision
region in terms of the responder counts $`y`$ that lead to practical
equivalence.

``` r

summary(des_onc)
```

The summary (not shown here) gives the first and last rows of the
calibration grid, along with the selected design point. These values can
be reported, e.g. as a table listing $`n^\ast`$, the ROPE region
$`\mathcal{R}_p = [0.18, 0.42]`$, the decision thresholds
$`\gamma_{\mathrm{eq}}, \gamma_{\mathrm{diff}}`$ and the resulting
ROPE-based power and type-I error. This is primarily helpful when
analyzing a specific design or the relationship of the operating
characteristics and the sample size.

### Visualization

We can inspect the operating characteristics and prior structure in more
detail.

``` r

plot(des_onc)
```

![Figure 8: Visualization of the calibrated ROPE-based oncology
single-arm one-stage phase II design with binary
endpoints.](figures/singlearm-onestage-rope-fig8.png)

Figure 8: Visualization of the calibrated ROPE-based oncology single-arm
one-stage phase II design with binary endpoints.

- The upper-left panel shows ROPE-based power and type-I error as
  functions of $`n`$.
- The upper-right panel summarizes the design numerically.
- The lower-left and middle panels overlay the design priors under
  $`H_0`$ and $`H_1`$.
- The lower-right panel shows the analysis prior.

For example, the equivalence design prior `Beta(36, 84)` reflects prior
belief that in realistic equivalence scenarios, the ORR is close to 30%,
whereas the non-equivalence design prior `Beta(60, 40)` reflects
scenarios with substantially higher ORR around 60%.

To see how the equivalence decision region changes with sample size, we
can plot the decision region directly:

``` r

plot(des_onc, what = "decision_region")
```

![Figure 9: Visualization of the equivalence region of ROPE-based
oncology single-arm one-stage phase II designs with binary endpoints for
increasing sample size.](figures/singlearm-onestage-rope-fig9.png)

Figure 9: Visualization of the equivalence region of ROPE-based oncology
single-arm one-stage phase II designs with binary endpoints for
increasing sample size.

This plot shows, for each evaluated sample size $`n`$, the range of
responder counts $`y`$ that would lead the trial to conclude practical
equivalence. For the selected $`n^\ast`$, this region is reported in the
upper right panel of Figure 8: If 20 to 35 patients show a success in
the phase II trial (out of $`n^\ast`$=94), then equivalence of the novel
drug or treatment to the reference probability $`p_0=0.30`$ (of the
standard of care) is established. Thus, we then accept
$`H_1:p \notin [\,p_0 - \delta,\; p_0 + \delta\,]`$.

Figure 8 also shows that both the Bayesian power and type-I-error rate
are calibrated.

## Example 2: Sensitivity analysis via grid exploration

Here we explore the impact of different design priors, ROPE half-widths
$`\delta`$ and the posterior probability threshold $`\gamma_{eq}`$ for
establishing equivalence.

    #> # A tibble: 9 × 6
    #>   delta gamma_eq n_star power_H1  type1_H0 feasible
    #>   <dbl>    <dbl>  <int>    <dbl>     <dbl> <lgl>   
    #> 1  0.1      0.75    138    0.818  0.000254 TRUE    
    #> 2  0.1      0.8     167    0.812  0.000111 TRUE    
    #> 3  0.1      0.9      NA   NA     NA        FALSE   
    #> 4  0.12     0.75     77    0.827  0.00200  TRUE    
    #> 5  0.12     0.8      94    0.823  0.000922 TRUE    
    #> 6  0.12     0.9     148    0.814  0.000156 TRUE    
    #> 7  0.15     0.75     41    0.817  0.0159   TRUE    
    #> 8  0.15     0.8      52    0.835  0.00769  TRUE    
    #> 9  0.15     0.9      78    0.820  0.00157  TRUE

| delta | gamma_eq | n_star | power_H1 | type1_H0 |
|------:|---------:|-------:|---------:|---------:|
|  0.10 |     0.75 |    138 |    0.818 |    0.000 |
|  0.10 |     0.80 |    167 |    0.812 |    0.000 |
|  0.10 |     0.90 |     NA |       NA |       NA |
|  0.12 |     0.75 |     77 |    0.827 |    0.002 |
|  0.12 |     0.80 |     94 |    0.823 |    0.001 |
|  0.12 |     0.90 |    148 |    0.814 |    0.000 |
|  0.15 |     0.75 |     41 |    0.817 |    0.016 |
|  0.15 |     0.80 |     52 |    0.835 |    0.008 |
|  0.15 |     0.90 |     78 |    0.820 |    0.002 |

Grid exploration for the oncology equivalence example: calibrated sample
size n\*, ROPE-based Bayesian power under H1, and ROPE-based Bayesian
type-I error under H0 for different ROPE half-widths and posterior
probability thresholds. {.table}

![Figure 10: Calibrated sample size n\* across ROPE widths and posterior
thresholds for the oncology equivalence phase II
trial.](figures/singlearm-onestage-rope-fig10.png)

Figure 10: Calibrated sample size n\* across ROPE widths and posterior
thresholds for the oncology equivalence phase II trial.

![Figure 10: ROPE-based Bayesian type-I-error at the calibrated sample
sizes for the oncology equivalence phase II trial for different ROPE
widths.](figures/singlearm-onestage-rope-fig11.png)

Figure 10: ROPE-based Bayesian type-I-error at the calibrated sample
sizes for the oncology equivalence phase II trial for different ROPE
widths.

------------------------------------------------------------------------

## Example 3: revisiting the first example with PCE(H0) and frequentist power

Here we revisit the first example of the oncology trial, now adding a
target constraint on the probability of compelling evidence for $`H_0`$
and also reporting frequentist power post-hoc for the resulting design:

``` r

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(knitr)

# Oncology-inspired equivalence example: revisited
n_min <- 10
n_max <- 300
p0    <- 0.30

# ROPE and evidence thresholds
delta      <- 0.12
gamma_eq   <- 0.925
gamma_diff <- 0.90

# Analysis prior
a <- 1
b <- 1

# Design priors as in the first example
da0 <- 60
db0 <- 40   # non-equivalence prior (H0)

da1 <- 36
db1 <- 84   # equivalence prior (H1)

# Calibration targets
target_power      <- 0.80   # Bayesian predictive power under H1
target_type1      <- 0.10   # Bayesian predictive type-I error under H0
target_pce_h0     <- 0.80   # predictive compelling evidence for H0
target_freq_power <- 0.80   # frequentist power at dp (here dp = p0)
target_freq_type1 <- 0.10   # frequentist type-I error at ROPE boundaries

# Point alternative for frequentist power
dp <- p0

# Design calibration in "full" mode
fit_pce_freq <- design_singlearm_onestage_rope(
  n_min      = n_min,
  n_max      = n_max,
  p0         = p0,
  delta      = delta,
  gamma_eq   = gamma_eq,
  gamma_diff = gamma_diff,
  direction  = "equivalence",
  a          = a,
  b          = b,
  da0        = da0,
  db0        = db0,
  da1        = da1,
  db1        = db1,
  calibration        = "full",
  dp                 = dp,
  target_power       = target_power,
  target_type1       = target_type1,
  target_pce_h0      = target_pce_h0,
  target_freq_power  = target_freq_power,
  target_freq_type1  = target_freq_type1,
  sustain_n          = 10,
  return_grid        = TRUE
)

fit_pce_freq
#> One-stage single-arm ROPE design
#> Direction: equivalence 
#> Calibration: full 
#> Search range n: 10 to 300 
#> Null probability p0: 0.3 
#> Margin delta: 0.12 
#> Probability threshold gamma_eq: 0.925 
#> Probability threshold gamma_diff: 0.9 
#> Analysis prior: Beta(1, 1)
#> Design prior (H0): Beta(60, 40)
#> Design prior (H1): Beta(36, 84)
#> Target Bayesian power: 0.8 
#> Target Bayesian type-I error: 0.1 
#> Target PCE(H0): 0.8 
#> Frequentist power point dp: 0.3 
#> Target frequentist power: 0.8 
#> Target frequentist type-I error: 0.1 
#> Sustain n: 10 
#> Selected sample size n*: 173 
#> Bayesian power(n*): 0.8166 
#> Bayesian type-I(n*): 0.0001 
#> PCE(H0)(n*): 0.9846 
#> Frequentist power(n*): 0.9597 
#> Frequentist type-I(n*): 0.0784 
#>  at p0 - delta: 0.0755 
#>  at p0 + delta: 0.0784 
#> Equivalence region: {39-63}
#> Compelling evidence for non-equivalence region: {0-24, 81-173}
```

You can summarise and visualise the calibrated design:

![Figure 12: Calibrated one-stage ROPE-based oncology equivalence phase
II design with additional constraints on the probability of compelling
evidence for the null hypothesis. In contrast to the earlier example,
the probability of compelling evidence must reach 80% now, and
frequentist power and type-I-error rate must also fulfill their
respective target constraints of 80% and
10%.](figures/singlearm-onestage-rope-fig12.png)

Figure 12: Calibrated one-stage ROPE-based oncology equivalence phase II
design with additional constraints on the probability of compelling
evidence for the null hypothesis. In contrast to the earlier example,
the probability of compelling evidence must reach 80% now, and
frequentist power and type-I-error rate must also fulfill their
respective target constraints of 80% and 10%.

``` r

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
# Extract selected row and key operating characteristics
sel <- fit_pce_freq$selected

summary_tab <- tibble(
  quantity = c(
    "Selected sample size n*",
    "Bayesian power under H1 at n*",
    "Bayesian type-I error under H0 at n*",
    "PCE(H0) at n*",
    "Frequentist power at p = p0",
    "Frequentist type-I error (worst boundary)"
  ),
  value = c(
    fit_pce_freq$n_star,
    sel$power,
    sel$type1,
    sel$pce_h0,
    sel$freq_power,
    sel$freq_type1
  )
)

kable(
  summary_tab,
  digits = 3,
  col.names = c("Quantity", "Value"),
  caption = "Operating characteristics of the calibrated equivalence design with constraints on Bayesian power, Bayesian type-I error, PCE(H0), and frequentist power/type-I error."
)
```

| Quantity                                  |   Value |
|:------------------------------------------|--------:|
| Selected sample size n\*                  | 173.000 |
| Bayesian power under H1 at n\*            |   0.817 |
| Bayesian type-I error under H0 at n\*     |   0.000 |
| PCE(H0) at n\*                            |   0.985 |
| Frequentist power at p = p0               |   0.960 |
| Frequentist type-I error (worst boundary) |   0.078 |

Operating characteristics of the calibrated equivalence design with
constraints on Bayesian power, Bayesian type-I error, PCE(H0), and
frequentist power/type-I error. {.table}

Optionally, you can compare this design to the original first example
(purely Bayesian calibration) by recomputing the first example and
putting both designs side by side in a small table:

``` r

des_onc_with_freq_power <- design_singlearm_onestage_rope(
  n_min = 20,
  n_max = 200,
  p0 = 0.30,
  delta = 0.12,
  gamma_eq = 0.80,
  
  # frequentist power at p = 0.3
  dp = 0.3,

  # Analysis prior p ~ Beta(a, b)
  a = 1, b = 1,

  # Design priors under H0 and H1
  da0 = 60, db0 = 40,   # H0: non-equivalence, mean ~0.60
  da1 = 36, db1 = 84,   # H1: equivalence, mean ~0.3

  target_power = 0.80,
  target_type1 = 0.10,
  target_freq_type1 = 0.10,
  target_freq_power = 0.80,
  sustain_n = 10,
  calibration = "Bayesian"
)

sel_orig <- des_onc_with_freq_power$selected
sel_new  <- fit_pce_freq$selected

comparison_tab <- tibble(
  design = c("Bayesian (original)", "Full (Bayes + frequentist + PCE(H0))"),
  n_star = c(sel_orig$n, fit_pce_freq$n),
  bayes_power = c(sel_orig$power, sel_new$power),
  bayes_type1 = c(sel_orig$type1, sel_new$type1),
  pce_h0 = c(sel_orig$pce_h0, sel_new$pce_h0),
  freq_power = c(sel_orig$freq_power, sel_new$freq_power),
  freq_type1 = c(sel_orig$freq_type1, sel_new$freq_type1)
)

kable(
  comparison_tab,
  digits = 3,
  caption = "Comparison of the original Bayesian calibration and the extended design with additional constraints on PCE(H0) and frequentist power/type-I error."
)
```

| design | n_star | bayes_power | bayes_type1 | pce_h0 | freq_power | freq_type1 |
|:---|---:|---:|---:|---:|---:|---:|
| Bayesian (original) | 94 | 0.823 | 0.001 | 0.973 | 0.925 | 0.240 |
| Full (Bayes + frequentist + PCE(H0)) | 173 | 0.817 | 0.000 | 0.985 | 0.960 | 0.078 |

Comparison of the original Bayesian calibration and the extended design
with additional constraints on PCE(H0) and frequentist power/type-I
error. {.table}

This third example stays within the equivalence framework but shows how
the **same posterior-threshold decision rule** can be calibrated to
satisfy additional Bayesian and frequentist criteria, including a lower
bound on predictive compelling evidence for $`H_0`$.

## Summary

This vignette has shown how
[`design_singlearm_onestage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_onestage_rope.md)
can be used to:

- define a baseline ROPE-based equivalence design in a realistic phase
  II range,

- quantify how the evidence threshold `gamma_eq`, ROPE width `delta`,
  design priors, and sustain requirement influence the required sample
  size, power, and type-I error.

In practice, we recommend exploring such grids of tuning parameters
collaboratively with clinicians, to arrive at a design where the ROPE
region, evidence thresholds, and priors are all clinically interpretable
and the resulting sample size is operationally feasible.

## References

Chen, Lichang, Jianhong Pan, Yanpeng Wu, et al. 2022. “Bayesian
Two-Stage Design for Phase II Oncology Trials with Binary Endpoint.”
*Statistics in Medicine* 41 (12): 2291–301.
<https://doi.org/10.1002/sim.9355>.

Hashim, Mahmoud, Talitha Vincken, Florint Kroi, et al. 2021. “A
Systematic Review of Noninferiority Margins in Oncology Clinical
Trials.” *Journal of Comparative Effectiveness Research* 10 (6): 443–55.
<https://doi.org/10.2217/cer-2020-0200>.

Hodges, J. L., and E. L. Lehmann. 1954. “Testing the Approximate
Validity of Statistical Hypotheses.” *Journal of the Royal Statistical
Society: Series B (Methodological)* 16 (2): 261–68.
<https://doi.org/10.1111/j.2517-6161.1954.tb00169.x>.

Kelter, Riko. 2021. “Bayesian Hodges-Lehmann Tests for Statistical
Equivalence in the Two-Sample Setting: Power Analysis, Type I Error
Rates and Equivalence Boundary Selection in Biomedical Research.” *BMC
Medical Research Methodology* 21 (171).
<https://doi.org/10.1186/s12874-021-01341-7>.

Kelter, Riko, and Alexander Schnurr. 2024. “The Bayesian
Group-Sequential Predictive Evidence Value Design for Phase II Clinical
Trials with Binary Endpoints.” *Statistics in Biosciences* 17: 442–78.
<https://doi.org/10.1007/s12561-024-09430-z>.

Kruschke, John K. 2014. *Doing Bayesian Data Analysis: A Tutorial with
R, JAGS, and Stan*. In *Doing Bayesian Data Analysis: A Tutorial with R,
JAGS, and Stan*, 2nd ed. Academic Press.
<https://doi.org/10.1016/B978-0-12-405888-0.09999-2>.

Kruschke, John K. 2018. “Rejecting or Accepting Parameter Values in
Bayesian Estimation.” *Advances in Methods and Practices in
Psychological Science* 1(2): 270–80.
<https://doi.org/10.1177/2515245918771304>.

Kruschke, John K., and T. M. Liddell. 2018. “The Bayesian New Statistics
: Hypothesis Testing, Estimation, Meta-Analysis, and Power Analysis from
a Bayesian Perspective.” *Psychonomic Bulletin and Review* 25: 178–206.
<https://doi.org/10.3758/s13423-016-1221-4>.

Lakens, Daniël, Anne M. Scheel, and Peder M. Isager. 2018. “Equivalence
Testing for Psychological Research: A Tutorial.” *Advances in Methods
and Practices in Psychological Science* 1 (2): 259–69.
<https://doi.org/10.1177/2515245918770963>.

Lee, Jack, and Diane D. Liu. 2008. “A Predictive Probability Design for
Phase II Cancer Clinical Trials.” *Clinical Trials* 5 (2): 93–106.
<https://doi.org/10.1177/1740774508089279>.

Liao, J. G., Vishal Midya, and Arthur Berg. 2020. “Connecting and
Contrasting the Bayes Factor and a Modified ROPE Procedure for Testing
Interval Null Hypotheses.” *American Statistician*, ahead of print.
<https://doi.org/10.1080/00031305.2019.1701550>.

Linde, Maximilian, Jorge N. Tendeiro, Ravi Selker, Eric Jan Wagenmakers,
and Don van Ravenzwaaij. 2023. “Decisions about Equivalence: A
Comparison of TOST, HDI-ROPE, and the Bayes Factor.” *Psychological
Methods* 28 (3): 740–55. <https://doi.org/10.1037/MET0000402>.

Pan, Jian, Rania Christoforou, Lucy Nives Wiedermann, and Marcel
Schweiker. 2025. “The Untapped Potential of Bayesian Region of Practical
Equivalence for Assessing Null Effects in Multi-Domain Research.”
*Building and Environment* 283: 113390.
<https://doi.org/10.1016/j.buildenv.2025.113390>.

Wellek, Stefan. 2010. *Testing Statistical Hypotheses of Equivalence and
Noninferiority*. In *Testing Statistical Hypotheses of Equivalence and
Noninferiority*. CRC Press. <https://doi.org/10.1201/ebk1439808184>.
