# Optimal full calibration for single-arm two-stage Bayes factor designs with binary endpoints

## Introduction: Background on hybrid and full calibration

In a Bayes-factor based design, decision rules are specified through
evidence thresholds for the Bayes factor $`BF_{01}`$ comparing a null
hypothesis $`H_0`$ to an alternative $`H_1`$. Calibration refers to how
these thresholds and the sample sizes are chosen so that the resulting
design has prespecified operating characteristics such as type‑I error
and power (Grieve 2022).

### Bayesian calibration

A **fully Bayesian calibration** chooses the design by controlling
Bayesian notions of type‑I error and power under explicit design priors
on the response probability. Both the false positive rate under $`H_0`$
and the power under $`H_1`$ are defined as prior‑predictive
probabilities with respect to these priors. This approach has been
developed for Bayes factors in binomial one- and two‑arm settings, where
closed‑form expressions and numerical integration can be used to obtain
sample sizes that satisfy Bayesian power and type‑I error targets
without Monte Carlo simulation, see (Kelter and Pawel 2025a), (Kelter
and Pawel 2025b) and (Kelter 2026). These methods provide a Bayesian
analogue of classical power analysis for Bayes factors and form the
basis of the one‑stage calibration routines in this package.

### Hybrid Bayes-frequentist calibration

At the same time, there is a long‑standing interest in **hybrid
Bayes–frequentist designs** that use Bayesian test statistics or
posterior quantities but demand that certain frequentist operating
characteristics are controlled (Grieve 2016). Examples include hybrid
designs in oncology trials (Lopez-Rey et al. 2025) and general
discussions about which metric to calibrate (Arjas and Gasbarra 2026)
(see also (Ryan et al. 2020)) in a Bayesian trial design. For a
systematic review about Bayesian sample size approaches in clinical
trials see (Marks et al. 2026), who found that the most common method
for sample size determination in Bayesian RCTs was a hybrid approach
(58% out of 164 clinical trials which explicitly used a Bayesian trial
design). Other examples include hybrid predictive monitoring schemes for
multi‑arm trials (Shi and Yin 2019), see also (Muehlemann et al. 2023).

Related ideas also appear in the literature on **calibrated Bayes
factors**, which seek Bayes factors whose behavior is aligned with
frequentist error control or prior‑predictive performance. This includes
multiplicity‑calibrated Bayesian hypothesis tests (Guo and Heitjan
2010), see also (Macrì Demartino et al. 2025).

### Current implementation

Against this backdrop, the **hybrid** and **full** calibration modes in
{bfbin2arm} implement Bayes–frequentist compromise designs for Bayes
factors in single‑ and two‑arm phase II trials with binary endpoints
(Kelter and Pawel 2025a, 2025b; Kelter 2026). In hybrid calibration,
power is defined in a Bayesian prior‑predictive sense under a design
prior for $`H_1`$, while type‑I error is controlled in a frequentist
sense at the null boundary $`p_0`$. In full calibration, both Bayesian
and frequentist error metrics are constrained simultaneously: Bayesian
power and Bayesian type‑I error under the design priors, and frequentist
power and type‑I error at fixed parameter values. This yields
Bayes-factor based designs that satisfy both sets of constraints and
make explicit the trade‑offs between Bayesian and frequentist notions of
error control in phase II trial design.

## Full Bayesian and frequentist calibration

Full calibration enforces both Bayesian and frequentist constraints
simultaneously for a single-arm two-stage design based on the Bayes
factor \\BF\_{01}\\.

In this mode, the design must satisfy:

- Bayesian power and Bayesian type-I error targets under specified
  design priors for \\H_1\\ and \\H_0\\,
- frequentist power and frequentist type-I error targets at a fixed
  point alternative \\dp\\ and at \\p = p_0\\.

This yields designs that simultaneously meet Bayesian planning criteria
and frequentist calibration requirements.

Full calibration is requested via

``` r

calibration = "full"
```

in
[`design_singlearm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_bf.md).
In this mode, feasibility requires:

- Bayesian power (prior-predictive under \\H_1\\) to be at least
  `target_power`,
- Bayesian type-I error (averaged under \\H_0\\) to be at most
  `target_type1`,
- frequentist power at \\p = dp\\ to be at least `target_freq_power`,
- frequentist type-I error at \\p = p_0\\ to be at most
  `target_freq_type1`.

The optional `target_ce_h0` can also be used to impose a lower bound on
the Bayesian probability of compelling evidence in favour of \\H_0\\,
computed under the \\H_0\\ design prior.

## Overview of the different calibration modes

The following table shows the different calibration modes available and
the use of the `power_cushion` parameter when calibrating the optimal
single-arm two-stage design:

[TABLE]

Use of power_cushion in the fixed-sample anchor step across calibration
modes. {.table .table .table-striped .table-condensed
style="font-size: 12px; "}

The first step of the calibration algorithm for an optimal single-arm
two-stage design therefore always consists in finding a sufficiently
large anchor sample size in the first step. Therefore, the relevant
power (frequentist or Bayesian) must exceed the target constraint plus
the power cushion specified as input. For Bayesian and full calibration,
the probability of compelling evidence (CE) must optionally also be
sufficiently large when searching the anchor sample size.

Based on this anchor sample size, step 2 of the algorithm consists in
computing the corrected operating characteristics for the two-stage
designs which are available based on the user input. The latter
includes, in particular, the minimum and maximum sample sizes `n1_min`
and `n2_max`, which limit the number of available single-arm two-stage
designs to analyze regarding optimality.

## Full calibration example

We first construct a fully calibrated design with moderate targets:

``` r

res_full <- design_singlearm_bf(
  n1_min = 5,
  n2_max = 100,
  k      = 1/10,
  k_f    = 3,
  p0     = 0.2,
  a0     = 1,
  b0     = 1,
  a1     = 1,
  b1     = 1,
  dp     = 0.5,
  da0    = 1,
  db0    = 1,
  da1    = 2.5,
  db1    = 2,
  type   = "direction",
  calibration       = "full",
  target_power      = 0.80,
  target_type1      = 0.05,
  target_freq_power = 0.80,
  target_freq_type1 = 0.1,
  power_cushion = 0.025
)
```

Note that we specified a power cushion via the parameter
`power_cushion = 0.025`. Otherwise, it might be impossible to find an
optimal design, for details on the underlying methodology see also the
vignette on hybrid calibration and the discussion section in (Kelter and
Pawel 2025b). We inspect the results produced by the calibration
algorithm:

``` r

summary(res_full)
#> Summary: Single-arm two-stage Bayes factor design
#> ---------------------------------------------------------
#> Feasible: TRUE
#> Calibration: full
#> Design prior under H0: Beta(1, 1) truncated to [0, p0]
#> Design prior under H1: Beta(2.5, 2) truncated to (p0, 1]
#> 
#> Selected design: n1 = 7, n2 = 24
#> 
#> Bayesian operating characteristics
#>   Power: 0.8072
#>   Type-I: 0.0043
#>   CE H0: NA
#>   EN H0: 9.88
#>   EN H1: 22.45
#> 
#> Frequentist operating characteristics
#>   Power: 0.8828
#>   Type-I: 0.0316
#>   EN H0: 14.20
#>   EN H1: 22.94
```

The summary reports:

- the selected interim and final sample sizes (`n1`, `n2`),
- Bayesian operating characteristics: power, type-I error, expected
  sample sizes, and (optionally) compelling evidence under \\H_0\\,
- frequentist operating characteristics: power and type-I error at
  \\dp\\ and \\p_0\\, expected sample sizes under \\H_0\\ and \\H_1\\.

A diagnostic plot illustrates how the search over interim sample sizes
balances these constraints:

``` r

if (isTRUE(res_full$feasible)) {
  plot(res_full)
}
```

![Figure 1: Output of the plot function for an optimal fully calibrated
single-arm two-stage design using Bayes factors. The top left panel
shows Bayesian and frequentist power, Bayesian type-I-error for varying
interim sample sizes. The top right panel provides information about the
optimal frequentist design found by the algorithm and its Bayesian and
frequentist operating characteristics. The lower left and right panels
visualize the analysis and design priors under the null and alternative
hypothesis. For the frequentist operating characteristics, these are
irrelevant. They influence only the Bayesian operating
characteristics.](figures/optimal_single_arm_two_stage_full_fig1.png)

Figure 1: Output of the plot function for an optimal fully calibrated
single-arm two-stage design using Bayes factors. The top left panel
shows Bayesian and frequentist power, Bayesian type-I-error for varying
interim sample sizes. The top right panel provides information about the
optimal frequentist design found by the algorithm and its Bayesian and
frequentist operating characteristics. The lower left and right panels
visualize the analysis and design priors under the null and alternative
hypothesis. For the frequentist operating characteristics, these are
irrelevant. They influence only the Bayesian operating characteristics.

## Why the design with $`n_1 = 8`$ is selected

In this example, several two-stage designs satisfy the full calibration
constraints. In particular, both the designs with $`n_1 = 6`$ and
$`n_1 = 8`$ meet the required Bayesian and frequentist power and type-I
error thresholds. Thus, the choice of the final design is not driven by
feasibility alone, but by the optimization criterion used after
feasibility has been established.

In the current implementation, the **optimal fully calibrated design**
is defined as the feasible design that minimizes the **Bayesian expected
sample size under $`H_0`$**. That is, among all designs satisfying the
Bayesian and frequentist operating-characteristic constraints, the
selected design is the one with the smallest value of `en_h0`.

For the present example, the design with $`n_1 = 6`$ is fully
calibrated, but its Bayesian expected sample size under $`H_0`$ is
``` math
\operatorname{EN}_{H_0}^{\mathrm{Bayes}} = 13.40.
```
The design with $`n_1 = 8`$ is also fully calibrated, but has the
smaller value
``` math
\operatorname{EN}_{H_0}^{\mathrm{Bayes}} = 11.09.
```
Because $`11.09 < 13.40`$, the design with $`n_1 = 8`$ is preferred and
is therefore returned as the optimal fully calibrated design (see also
the upper right panel in Figure 1 which reports that expected sample
size of the optimal design).

This behavior reflects the current philosophy of the package:
calibration constraints determine which designs are admissible, and
among these admissible designs the optimization is performed with
respect to the Bayesian expected sample size under $`H_0`$. Intuitively,
this favors designs that stop early more efficiently when the null
hypothesis is true, while still maintaining the required Bayesian and
frequentist error control.

Possible future extensions could allow the optimization criterion to be
changed. For example, one might instead minimize the **frequentist**
expected sample size under $`H_0`$, or define a compromise criterion
based on a weighted average of the Bayesian and frequentist expected
sample sizes under $`H_0`$, such as
``` math
w \cdot \operatorname{EN}_{H_0}^{\mathrm{Bayes}}
+
(1-w) \cdot \operatorname{EN}_{H_0}^{\mathrm{Freq}},
\qquad 0 \le w \le 1.
```
Such extensions would make it possible to tailor the notion of
optimality more closely to the user’s preferred balance between Bayesian
and frequentist design perspectives.

## Comparison of calibration modes

It is often instructive to compare the designs obtained under different
calibration modes with identical thresholds and priors.

``` r

res_bayes <- design_singlearm_bf(
  n1_min = 5,
  n2_max = 100,
  k      = 1/10,
  k_f    = 3,
  p0     = 0.2,
  a0     = 1,
  b0     = 1,
  a1     = 1,
  b1     = 1,
  dp     = 0.5,
  da0    = 1,
  db0    = 1,
  da1    = 2.5,
  db1    = 2,
  type   = "direction",
  calibration       = "Bayesian",
  target_power      = 0.80,
  target_type1      = 0.05,
  target_freq_power = 0.80,
  target_freq_type1 = 0.05,
  power_cushion = 0.025
)

res_freq <- design_singlearm_bf(
  n1_min = 5,
  n2_max = 100,
  k      = 1/10,
  k_f    = 3,
  p0     = 0.2,
  a0     = 1,
  b0     = 1,
  a1     = 1,
  b1     = 1,
  dp     = 0.5,
  da0    = 1,
  db0    = 1,
  da1    = 2.5,
  db1    = 2,
  type   = "direction",
  calibration       = "frequentist",
  target_power      = 0.80,
  target_type1      = 0.05,
  target_freq_power = 0.80,
  target_freq_type1 = 0.05,
  power_cushion = 0.025
)

res_hybrid <- design_singlearm_bf(
  n1_min = 5,
  n2_max = 100,
  k      = 1/10,
  k_f    = 3,
  p0     = 0.2,
  a0     = 1,
  b0     = 1,
  a1     = 1,
  b1     = 1,
  dp     = 0.5,
  da0    = 1,
  db0    = 1,
  da1    = 2.5,
  db1    = 2,
  type   = "direction",
  calibration       = "hybrid",
  target_power      = 0.80,
  target_type1      = 0.05,
  target_freq_power = 0.80,
  target_freq_type1 = 0.05,
  power_cushion = 0.025
)
```

We then compare the main operating characteristics:

|  | calibration | n1 | n2 | bayes_power | bayes_type1 | freq_power | freq_type1 | bayes_en_h0 | bayes_en_h1 | freq_en_h0 | freq_en_h1 |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| n1 | Bayesian | 7 | 24 | 0.807 | 0.004 | 0.883 | 0.032 | 9.881 | 22.451 | 14.196 | 22.938 |
| n11 | frequentist | 7 | 17 | 0.775 | 0.006 | 0.812 | 0.035 | 8.695 | 16.089 | 11.233 | 16.375 |
| n12 | hybrid | 7 | 24 | 0.807 | 0.004 | 0.883 | 0.032 | 9.881 | 22.451 | 14.196 | 22.938 |
| n13 | full | 7 | 24 | 0.807 | 0.004 | 0.883 | 0.032 | 9.881 | 22.451 | 14.196 | 22.938 |

This comparison table highlights:

- how the chosen calibration mode influences `n1`, `n2`,
- the trade-offs between Bayesian and frequentist power and type-I
  error,
- differences in expected sample size under \\H_0\\ and \\H_1\\.

The full calibration design will usually be more demanding in terms of
sample size than the purely Bayesian, purely frequentist, or hybrid
designs, because it must simultaneously satisfy all constraints. In this
specific example, the Bayesian and hybrid optimal single-arm two-stage
design is fully calibrated. That is, it also satisfies the frequentist
target constraints. We can also check this manually as follows:

``` r

summary(res_bayes)
#> Summary: Single-arm two-stage Bayes factor design
#> ---------------------------------------------------------
#> Feasible: TRUE
#> Calibration: Bayesian
#> Design prior under H0: Beta(1, 1) truncated to [0, p0]
#> Design prior under H1: Beta(2.5, 2) truncated to (p0, 1]
#> 
#> Selected design: n1 = 7, n2 = 24
#> 
#> Bayesian operating characteristics
#>   Power: 0.8072
#>   Type-I: 0.0043
#>   CE H0: NA
#>   EN H0: 9.88
#>   EN H1: 22.45
#> 
#> Frequentist operating characteristics
#>   Power: 0.8828
#>   Type-I: 0.0316
#>   EN H0: 14.20
#>   EN H1: 22.94
```

which shows that the frequentist power and type-I-error rate meet our
target constraints. Likewise, we see based on

``` r

summary(res_hybrid)
#> Summary: Single-arm two-stage Bayes factor design
#> ---------------------------------------------------------
#> Feasible: TRUE
#> Calibration: hybrid
#> Design prior under H0: Beta(1, 1) truncated to [0, p0]
#> Design prior under H1: Beta(2.5, 2) truncated to (p0, 1]
#> 
#> Selected design: n1 = 7, n2 = 24
#> 
#> Bayesian operating characteristics
#>   Power: 0.8072
#>   Type-I: 0.0043
#>   CE H0: NA
#>   EN H0: 9.88
#>   EN H1: 22.45
#> 
#> Frequentist operating characteristics
#>   Power: 0.8828
#>   Type-I: 0.0316
#>   EN H0: 14.20
#>   EN H1: 22.94
```

that the same holds for the optimal hybrid single-arm two-stage design
(the latter is identical to the Bayesian one in this case).

In contrast, the frequentist optimal single-arm two-stage design has
different sample sizes $`n_1`$ and $`n_2`$, and also different operating
characteristics:

``` r

summary(res_freq)
#> Summary: Single-arm two-stage Bayes factor design
#> ---------------------------------------------------------
#> Feasible: TRUE
#> Calibration: frequentist
#> Design prior under H0: Beta(1, 1) truncated to [0, p0]
#> Design prior under H1: Beta(2.5, 2) truncated to (p0, 1]
#> 
#> Selected design: n1 = 7, n2 = 17
#> 
#> Bayesian operating characteristics
#>   Power: 0.7752
#>   Type-I: 0.0056
#>   CE H0: NA
#>   EN H0: 8.69
#>   EN H1: 16.09
#> 
#> Frequentist operating characteristics
#>   Power: 0.8119
#>   Type-I: 0.0351
#>   EN H0: 11.23
#>   EN H1: 16.38
```

## Practical recommendations for full calibration

Full calibration is most appropriate when:

- both Bayesian and frequentist perspectives must be satisfied, for
  example, to reconcile Bayesian planning with regulatory frequentist
  criteria, and
- there is sufficient flexibility in sample size (i.e. reasonably large
  `n2_max`) so that a design can exist that meets all targets.

If no feasible design is found, consider:

- relaxing at least one of the targets (often the Bayesian or
  frequentist power),
- increasing `n2_max` to allow more patient numbers,
- or reverting to a hybrid calibration if Bayesian power is most
  important but a strict frequentist power target at a single point
  alternative is not essential.

## References

Arjas, Elja, and Dario Gasbarra. 2026. *Is Control of Type I Error Rate
Needed in Bayesian Clinical Trial Designs?* arXiv:2312.15222. arXiv.
<https://doi.org/10.48550/arXiv.2312.15222>.

Grieve, Andrew P. 2016. “Idle Thoughts of a ’Well-Calibrated’ Bayesian
in Clinical Drug Development.” *Pharmaceutical Statistics* 15 (2):
96–108. <https://doi.org/10.1002/PST.1736>.

Grieve, Andrew P. 2022. *Hybrid Frequentist/Bayesian Power and Bayesian
Power in Planning and Clinical Trials*. Chapman & Hall, CRC Press.

Guo, Mengye, and Daniel F. Heitjan. 2010. “Multiplicity-Calibrated
Bayesian Hypothesis Tests.” *Biostatistics (Oxford, England)* 11 (3):
473–83. <https://doi.org/10.1093/biostatistics/kxq012>.

Kelter, Riko. 2026. *Power and Sample Size Calculations for Bayes
Factors in Two-Arm Clinical Phase II Trials with Binary Endpoints*.
<https://arxiv.org/abs/2603.01715>.

Kelter, Riko, and Samuel Pawel. 2025a. *Bayesian Power and Sample Size
Calculations for Bayes Factors in the Binomial Setting*.
<https://arxiv.org/abs/2502.02914>.

Kelter, Riko, and Samuel Pawel. 2025b. *The Bayesian Optimal Two-Stage
Design for Clinical Phase II Trials Based on Bayes Factors*.
<https://arxiv.org/abs/2511.23144>.

Lopez-Rey, Borja G., Gerard Carot-Sans, Dan Ouchi, Ferran Torres, and
Caridad Pontes. 2025. “Use of Bayesian Approaches in Oncology Clinical
Trials: A Cross-Sectional Analysis.” *Frontiers in Pharmacology* 16
(March). <https://doi.org/10.3389/fphar.2025.1548997>.

Macrì Demartino, Roberto, Leonardo Egidi, Nicola Torelli, and Ioannis
Ntzoufras. 2025. “Eliciting Prior Information from Clinical Trials via
Calibrated Bayes Factor.” *Computational Statistics & Data Analysis* 209
(September): 108180. <https://doi.org/10.1016/j.csda.2025.108180>.

Marks, Yanara, Jessie Cunningham, Arlene Jiang, et al. 2026. “A
Systematic Review of Sample Size Determination in Bayesian Randomized
Clinical Trials: Full Bayesian Methods Are Rarely Used.” *BMC Medical
Research Methodology*, ahead of print, April.
<https://doi.org/10.1186/s12874-026-02854-9>.

Muehlemann, Natalia, Tianjian Zhou, Rajat Mukherjee, Munshi Imran
Hossain, Satrajit Roychoudhury, and Estelle Russek-Cohen. 2023. “A
Tutorial on Modern Bayesian Methods in Clinical Trials.” *Therapeutic
Innovation & Regulatory Science* 57 (3): 402–16.
<https://doi.org/10.1007/s43441-023-00515-3>.

Ryan, Elizabeth G., Kristian Brock, Simon Gates, and Daniel Slade. 2020.
“Do We Need to Adjust for Interim Analyses in a Bayesian Adaptive Trial
Design?” *BMC Medical Research Methodology* 20 (1).
<https://doi.org/10.1186/S12874-020-01042-7>.

Shi, Haolun, and Guosheng Yin. 2019. “Control of Type I Error Rates in
Bayesian Sequential Designs.” *Bayesian Analysis* 14 (2): 399–425.
<https://doi.org/10.1214/18-BA1109>.
