# Sample size calibration for two-arm binomial Bayes factor designs

Searches over a grid of total sample sizes n to find the smallest n such
that Bayesian power, Bayesian type-I error, and probability of
compelling evidence under H0 meet specified design criteria. Optionally,
frequentist type-I error and power constraints are also evaluated.
Unequal fixed randomisation between the two arms is allowed via alloc1
and alloc2.

Backward-compatible wrapper around
[`design_twoarm_onestage_bf()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_bf.md).

## Usage

``` r
ntwoarmbinbf01(
  k = 1/3,
  k_f = 3,
  power = 0.8,
  alpha = 0.05,
  pce_H0 = 0.9,
  test = c("BF01", "BF+0", "BF-0", "BF+-"),
  nrange = c(10, 150),
  n_step = 1,
  progress = TRUE,
  compute_freq_t1e = FALSE,
  p1_grid = seq(0.01, 0.99, 0.02),
  p2_grid = seq(0.01, 0.99, 0.02),
  p1_power = NULL,
  p2_power = NULL,
  a_0_d = 1,
  b_0_d = 1,
  a_0_a = 1,
  b_0_a = 1,
  a_1_d = 1,
  b_1_d = 1,
  a_2_d = 1,
  b_2_d = 1,
  a_1_a = 1,
  b_1_a = 1,
  a_2_a = 1,
  b_2_a = 1,
  output = c("plot", "numeric"),
  a_1_d_Hminus = 1,
  b_1_d_Hminus = 1,
  a_2_d_Hminus = 1,
  b_2_d_Hminus = 1,
  a_1_a_Hminus = 1,
  b_1_a_Hminus = 1,
  a_2_a_Hminus = 1,
  b_2_a_Hminus = 1,
  alloc1 = 0.5,
  alloc2 = 0.5,
  sustain_n = 10L
)

ntwoarmbinbf01(
  k = 1/3,
  k_f = 3,
  power = 0.8,
  alpha = 0.05,
  pce_H0 = 0.9,
  test = c("BF01", "BF+0", "BF-0", "BF+-"),
  nrange = c(10, 150),
  n_step = 1,
  progress = TRUE,
  compute_freq_t1e = FALSE,
  p1_grid = seq(0.01, 0.99, 0.02),
  p2_grid = seq(0.01, 0.99, 0.02),
  p1_power = NULL,
  p2_power = NULL,
  a_0_d = 1,
  b_0_d = 1,
  a_0_a = 1,
  b_0_a = 1,
  a_1_d = 1,
  b_1_d = 1,
  a_2_d = 1,
  b_2_d = 1,
  a_1_a = 1,
  b_1_a = 1,
  a_2_a = 1,
  b_2_a = 1,
  output = c("plot", "numeric"),
  a_1_d_Hminus = 1,
  b_1_d_Hminus = 1,
  a_2_d_Hminus = 1,
  b_2_d_Hminus = 1,
  a_1_a_Hminus = 1,
  b_1_a_Hminus = 1,
  a_2_a_Hminus = 1,
  b_2_a_Hminus = 1,
  alloc1 = 0.5,
  alloc2 = 0.5,
  sustain_n = 10L
)
```

## Arguments

- k:

  Evidence threshold for rejecting the null (inverted BF).

- k_f:

  Evidence threshold for "compelling evidence" in favour of the null.

- power:

  Target Bayesian power.

- alpha:

  Target Bayesian type-I error.

- pce_H0:

  Target probability of compelling evidence under \\H_0\\.

- test:

  Character string, one of `"BF01"`, `"BF+0"`, `"BF-0"`, `"BF+-"`.

- nrange:

  Integer vector of length 2 giving the search range for total n.

- n_step:

  Step size for n. Currently only `n_step = 1` is supported in the
  object-based calibration workflow.

- progress:

  Logical; if `TRUE`, print progress to the console.

- compute_freq_t1e:

  Logical; if `TRUE`, compute frequentist type-I error over a grid.

- p1_grid, p2_grid:

  Grids of true proportions for frequentist T1E.

- p1_power, p2_power:

  Optional true proportions for frequentist power.

- a_0_d, b_0_d, a_0_a, b_0_a:

  Shape parameters for design and analysis priors under \\H_0\\.

- a_1_d, b_1_d, a_2_d, b_2_d:

  Shape parameters for design priors under \\H_1\\ or \\H\_+\\.

- a_1_a, b_1_a, a_2_a, b_2_a:

  Shape parameters for analysis priors under \\H_1\\ or \\H\_+\\.

- output:

  `"plot"` or `"numeric"`.

- a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus:

  Optional design priors under \\H\_-\\ for directional tests.

- a_1_a_Hminus, b_1_a_Hminus, a_2_a_Hminus, b_2_a_Hminus:

  Shape parameters for analysis priors under H-.

- alloc1, alloc2:

  Fixed randomisation probabilities for arm 1 and arm 2; must be
  positive and sum to 1.

- sustain_n:

  Non-negative integer. A candidate total sample size is considered
  feasible only if the relevant target constraints hold at that total
  sample size and for the next `sustain_n` larger total sample sizes in
  the search range.

## Value

If `output = "plot"`, returns invisibly a list with recommended sample
sizes and a ggplot object printed to the device. If
`output = "numeric"`, returns a list with recommended n and summary.

If `output = "numeric"`, returns a `"twoarm_onestage_bf_design"` object.
If `output = "plot"`, the plot is printed and the design object is
returned invisibly.

## Examples

``` r
# Standard calibration with equal allocation: power 80%, type-I 5%, CE(H0) 80%
# \donttest{
ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.8, output = "numeric")
#> n_total=10 (n1=5, n2=5): power=0.333, type1=0.058, ce=0.182
#> n_total=11 (n1=6, n2=5): power=0.286, type1=0.034, ce=0.167
#> n_total=12 (n1=6, n2=6): power=0.245, type1=0.018, ce=0.154
#> n_total=13 (n1=6, n2=7): power=0.357, type1=0.049, ce=0.143
#> n_total=14 (n1=7, n2=7): power=0.312, type1=0.028, ce=0.133
#> n_total=15 (n1=8, n2=7): power=0.306, type1=0.023, ce=0.125
#> n_total=16 (n1=8, n2=8): power=0.370, type1=0.039, ce=0.118
#> n_total=17 (n1=8, n2=9): power=0.356, type1=0.030, ce=0.111
#> n_total=18 (n1=9, n2=9): power=0.380, type1=0.035, ce=0.105
#> n_total=19 (n1=10, n2=9): power=0.400, type1=0.037, ce=0.100
#> n_total=20 (n1=10, n2=10): power=0.380, type1=0.027, ce=0.145
#> n_total=21 (n1=10, n2=11): power=0.409, type1=0.035, ce=0.186
#> n_total=22 (n1=11, n2=11): power=0.417, type1=0.033, ce=0.219
#> n_total=23 (n1=12, n2=11): power=0.397, type1=0.025, ce=0.210
#> n_total=24 (n1=12, n2=12): power=0.450, type1=0.039, ce=0.202
#> n_total=25 (n1=12, n2=13): power=0.429, type1=0.029, ce=0.194
#> n_total=26 (n1=13, n2=13): power=0.429, type1=0.027, ce=0.187
#> n_total=27 (n1=14, n2=13): power=0.457, type1=0.033, ce=0.238
#> n_total=28 (n1=14, n2=14): power=0.436, type1=0.024, ce=0.255
#> n_total=29 (n1=14, n2=15): power=0.467, type1=0.032, ce=0.221
#> n_total=30 (n1=15, n2=15): power=0.461, type1=0.027, ce=0.238
#> n_total=31 (n1=16, n2=15): power=0.463, type1=0.027, ce=0.231
#> n_total=32 (n1=16, n2=16): power=0.471, type1=0.026, ce=0.245
#> n_total=33 (n1=16, n2=17): power=0.471, type1=0.024, ce=0.257
#> n_total=34 (n1=17, n2=17): power=0.494, type1=0.030, ce=0.268
#> n_total=35 (n1=18, n2=17): power=0.491, type1=0.027, ce=0.260
#> n_total=36 (n1=18, n2=18): power=0.488, type1=0.025, ce=0.270
#> n_total=37 (n1=18, n2=19): power=0.505, type1=0.028, ce=0.279
#> n_total=38 (n1=19, n2=19): power=0.495, type1=0.023, ce=0.286
#> n_total=39 (n1=20, n2=19): power=0.519, type1=0.030, ce=0.279
#> n_total=40 (n1=20, n2=20): power=0.512, type1=0.026, ce=0.286
#> n_total=41 (n1=20, n2=21): power=0.515, type1=0.025, ce=0.292
#> n_total=42 (n1=21, n2=21): power=0.529, type1=0.028, ce=0.297
#> n_total=43 (n1=22, n2=21): power=0.522, type1=0.024, ce=0.326
#> n_total=44 (n1=22, n2=22): power=0.537, type1=0.028, ce=0.318
#> n_total=45 (n1=22, n2=23): power=0.533, type1=0.025, ce=0.355
#> n_total=46 (n1=23, n2=23): power=0.535, type1=0.025, ce=0.336
#> n_total=47 (n1=24, n2=23): power=0.543, type1=0.026, ce=0.349
#> n_total=48 (n1=24, n2=24): power=0.541, type1=0.024, ce=0.390
#> n_total=49 (n1=24, n2=25): power=0.554, type1=0.027, ce=0.410
#> n_total=50 (n1=25, n2=25): power=0.553, type1=0.026, ce=0.428
#> n_total=51 (n1=26, n2=25): power=0.547, type1=0.022, ce=0.436
#> n_total=52 (n1=26, n2=26): power=0.554, type1=0.023, ce=0.448
#> n_total=53 (n1=26, n2=27): power=0.553, type1=0.022, ce=0.491
#> n_total=54 (n1=27, n2=27): power=0.566, type1=0.025, ce=0.481
#> n_total=55 (n1=28, n2=27): power=0.559, type1=0.022, ce=0.502
#> n_total=56 (n1=28, n2=28): power=0.568, type1=0.024, ce=0.524
#> n_total=57 (n1=28, n2=29): power=0.568, type1=0.022, ce=0.511
#> n_total=58 (n1=29, n2=29): power=0.567, type1=0.021, ce=0.528
#> n_total=59 (n1=30, n2=29): power=0.576, type1=0.023, ce=0.525
#> n_total=60 (n1=30, n2=30): power=0.572, type1=0.021, ce=0.532
#> n_total=61 (n1=30, n2=31): power=0.583, type1=0.023, ce=0.561
#> n_total=62 (n1=31, n2=31): power=0.578, type1=0.021, ce=0.546
#> n_total=63 (n1=32, n2=31): power=0.581, type1=0.021, ce=0.559
#> n_total=64 (n1=32, n2=32): power=0.588, type1=0.022, ce=0.547
#> n_total=65 (n1=32, n2=33): power=0.583, type1=0.020, ce=0.557
#> n_total=66 (n1=33, n2=33): power=0.593, type1=0.022, ce=0.569
#> n_total=67 (n1=34, n2=33): power=0.590, type1=0.020, ce=0.565
#> n_total=68 (n1=34, n2=34): power=0.599, type1=0.022, ce=0.583
#> n_total=69 (n1=34, n2=35): power=0.597, type1=0.021, ce=0.582
#> n_total=70 (n1=35, n2=35): power=0.597, type1=0.020, ce=0.577
#> n_total=71 (n1=36, n2=35): power=0.604, type1=0.021, ce=0.592
#> n_total=72 (n1=36, n2=36): power=0.599, type1=0.019, ce=0.607
#> n_total=73 (n1=36, n2=37): power=0.609, type1=0.021, ce=0.609
#> n_total=74 (n1=37, n2=37): power=0.604, type1=0.019, ce=0.617
#> n_total=75 (n1=38, n2=37): power=0.611, type1=0.021, ce=0.645
#> n_total=76 (n1=38, n2=38): power=0.611, type1=0.020, ce=0.619
#> n_total=77 (n1=38, n2=39): power=0.612, type1=0.019, ce=0.647
#> n_total=78 (n1=39, n2=39): power=0.616, type1=0.020, ce=0.635
#> n_total=79 (n1=40, n2=39): power=0.615, type1=0.019, ce=0.645
#> n_total=80 (n1=40, n2=40): power=0.623, type1=0.021, ce=0.649
#> n_total=81 (n1=40, n2=41): power=0.619, type1=0.019, ce=0.650
#> n_total=82 (n1=41, n2=41): power=0.628, type1=0.021, ce=0.669
#> n_total=83 (n1=42, n2=41): power=0.626, type1=0.020, ce=0.648
#> n_total=84 (n1=42, n2=42): power=0.623, type1=0.018, ce=0.674
#> n_total=85 (n1=42, n2=43): power=0.629, type1=0.020, ce=0.649
#> n_total=86 (n1=43, n2=43): power=0.627, type1=0.018, ce=0.669
#> n_total=87 (n1=44, n2=43): power=0.634, type1=0.020, ce=0.649
#> n_total=88 (n1=44, n2=44): power=0.631, type1=0.018, ce=0.665
#> n_total=89 (n1=44, n2=45): power=0.637, type1=0.020, ce=0.652
#> n_total=90 (n1=45, n2=45): power=0.637, type1=0.019, ce=0.660
#> n_total=91 (n1=46, n2=45): power=0.636, type1=0.018, ce=0.669
#> n_total=92 (n1=46, n2=46): power=0.639, type1=0.019, ce=0.655
#> n_total=93 (n1=46, n2=47): power=0.638, type1=0.018, ce=0.697
#> n_total=94 (n1=47, n2=47): power=0.643, type1=0.019, ce=0.651
#> n_total=95 (n1=48, n2=47): power=0.642, type1=0.018, ce=0.698
#> n_total=96 (n1=48, n2=48): power=0.649, type1=0.020, ce=0.662
#> n_total=97 (n1=48, n2=49): power=0.645, type1=0.018, ce=0.694
#> n_total=98 (n1=49, n2=49): power=0.644, type1=0.017, ce=0.691
#> n_total=99 (n1=50, n2=49): power=0.648, type1=0.018, ce=0.693
#> n_total=100 (n1=50, n2=50): power=0.647, type1=0.017, ce=0.714
#> n_total=101 (n1=50, n2=51): power=0.652, type1=0.018, ce=0.698
#> n_total=102 (n1=51, n2=51): power=0.651, type1=0.017, ce=0.731
#> n_total=103 (n1=52, n2=51): power=0.656, type1=0.018, ce=0.702
#> n_total=104 (n1=52, n2=52): power=0.651, type1=0.016, ce=0.731
#> n_total=105 (n1=52, n2=53): power=0.653, type1=0.016, ce=0.706
#> n_total=106 (n1=53, n2=53): power=0.656, type1=0.017, ce=0.727
#> n_total=107 (n1=54, n2=53): power=0.657, type1=0.017, ce=0.710
#> n_total=108 (n1=54, n2=54): power=0.659, type1=0.017, ce=0.724
#> n_total=109 (n1=54, n2=55): power=0.660, type1=0.017, ce=0.714
#> n_total=110 (n1=55, n2=55): power=0.663, type1=0.017, ce=0.720
#> n_total=111 (n1=56, n2=55): power=0.662, type1=0.016, ce=0.719
#> n_total=112 (n1=56, n2=56): power=0.665, type1=0.017, ce=0.720
#> n_total=113 (n1=56, n2=57): power=0.664, type1=0.016, ce=0.731
#> n_total=114 (n1=57, n2=57): power=0.664, type1=0.016, ce=0.716
#> n_total=115 (n1=58, n2=57): power=0.668, type1=0.017, ce=0.740
#> n_total=116 (n1=58, n2=58): power=0.667, type1=0.016, ce=0.713
#> n_total=117 (n1=58, n2=59): power=0.671, type1=0.017, ce=0.742
#> n_total=118 (n1=59, n2=59): power=0.670, type1=0.016, ce=0.709
#> n_total=119 (n1=60, n2=59): power=0.674, type1=0.017, ce=0.740
#> n_total=120 (n1=60, n2=60): power=0.672, type1=0.016, ce=0.706
#> n_total=121 (n1=60, n2=61): power=0.672, type1=0.015, ce=0.739
#> n_total=122 (n1=61, n2=61): power=0.675, type1=0.016, ce=0.749
#> n_total=123 (n1=62, n2=61): power=0.674, type1=0.015, ce=0.736
#> n_total=124 (n1=62, n2=62): power=0.679, type1=0.016, ce=0.761
#> n_total=125 (n1=62, n2=63): power=0.677, type1=0.015, ce=0.734
#> n_total=126 (n1=63, n2=63): power=0.682, type1=0.017, ce=0.761
#> n_total=127 (n1=64, n2=63): power=0.679, type1=0.015, ce=0.732
#> n_total=128 (n1=64, n2=64): power=0.682, type1=0.016, ce=0.761
#> n_total=129 (n1=64, n2=65): power=0.682, type1=0.015, ce=0.745
#> n_total=130 (n1=65, n2=65): power=0.681, type1=0.015, ce=0.761
#> n_total=131 (n1=66, n2=65): power=0.684, type1=0.016, ce=0.750
#> n_total=132 (n1=66, n2=66): power=0.683, type1=0.015, ce=0.761
#> n_total=133 (n1=66, n2=67): power=0.687, type1=0.016, ce=0.756
#> n_total=134 (n1=67, n2=67): power=0.686, type1=0.015, ce=0.760
#> n_total=135 (n1=68, n2=67): power=0.690, type1=0.016, ce=0.762
#> n_total=136 (n1=68, n2=68): power=0.688, type1=0.015, ce=0.757
#> n_total=137 (n1=68, n2=69): power=0.689, type1=0.015, ce=0.769
#> n_total=138 (n1=69, n2=69): power=0.690, type1=0.015, ce=0.755
#> n_total=139 (n1=70, n2=69): power=0.690, type1=0.015, ce=0.771
#> n_total=140 (n1=70, n2=70): power=0.692, type1=0.015, ce=0.752
#> n_total=141 (n1=70, n2=71): power=0.692, type1=0.015, ce=0.772
#> n_total=142 (n1=71, n2=71): power=0.695, type1=0.015, ce=0.751
#> n_total=143 (n1=72, n2=71): power=0.694, type1=0.015, ce=0.772
#> n_total=144 (n1=72, n2=72): power=0.698, type1=0.016, ce=0.748
#> n_total=145 (n1=72, n2=73): power=0.696, type1=0.015, ce=0.771
#> n_total=146 (n1=73, n2=73): power=0.697, type1=0.015, ce=0.774
#> n_total=147 (n1=74, n2=73): power=0.698, type1=0.015, ce=0.771
#> n_total=148 (n1=74, n2=74): power=0.698, type1=0.015, ce=0.783
#> n_total=149 (n1=74, n2=75): power=0.700, type1=0.015, ce=0.770
#> n_total=150 (n1=75, n2=75): power=0.700, type1=0.015, ce=0.785
#> 
#> One-stage two-arm Bayes factor design
#> ------------------------------------
#> Mode: optimal
#> Status: No feasible one-stage two-arm design found.
#> Calibration: Bayesian
#> Optional freq. Type-I reporting: off
#> Design: no feasible design found
# }

# 1:2 allocation (control:treatment) via alloc1 = 1/3, alloc2 = 2/3
# \donttest{
ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.8,
               alloc1 = 1/3, alloc2 = 2/3, output = "numeric")
#> n_total=10 (n1=3, n2=7): power=0.250, type1=0.035, ce=0.000
#> n_total=11 (n1=4, n2=7): power=0.250, type1=0.029, ce=0.167
#> n_total=12 (n1=4, n2=8): power=0.311, type1=0.044, ce=0.154
#> n_total=13 (n1=4, n2=9): power=0.280, type1=0.031, ce=0.143
#> n_total=14 (n1=5, n2=9): power=0.300, type1=0.029, ce=0.133
#> n_total=15 (n1=5, n2=10): power=0.303, type1=0.029, ce=0.125
#> n_total=16 (n1=5, n2=11): power=0.333, type1=0.036, ce=0.118
#> n_total=17 (n1=6, n2=11): power=0.333, type1=0.029, ce=0.183
#> n_total=18 (n1=6, n2=12): power=0.352, type1=0.032, ce=0.175
#> n_total=19 (n1=6, n2=13): power=0.367, type1=0.036, ce=0.168
#> n_total=20 (n1=7, n2=13): power=0.375, type1=0.031, ce=0.157
#> n_total=21 (n1=7, n2=14): power=0.383, type1=0.033, ce=0.152
#> n_total=22 (n1=7, n2=15): power=0.391, type1=0.034, ce=0.146
#> n_total=23 (n1=8, n2=15): power=0.403, type1=0.032, ce=0.138
#> n_total=24 (n1=8, n2=16): power=0.405, type1=0.031, ce=0.133
#> n_total=25 (n1=8, n2=17): power=0.395, type1=0.026, ce=0.129
#> n_total=26 (n1=9, n2=17): power=0.422, type1=0.030, ce=0.192
#> n_total=27 (n1=9, n2=18): power=0.421, type1=0.028, ce=0.186
#> n_total=28 (n1=9, n2=19): power=0.440, type1=0.033, ce=0.211
#> n_total=29 (n1=10, n2=19): power=0.436, type1=0.028, ce=0.201
#> n_total=30 (n1=10, n2=20): power=0.450, type1=0.031, ce=0.222
#> n_total=31 (n1=10, n2=21): power=0.446, type1=0.029, ce=0.217
#> n_total=32 (n1=11, n2=21): power=0.462, type1=0.030, ce=0.207
#> n_total=33 (n1=11, n2=22): power=0.464, type1=0.029, ce=0.202
#> n_total=34 (n1=11, n2=23): power=0.465, type1=0.029, ce=0.198
#> n_total=35 (n1=12, n2=23): power=0.481, type1=0.031, ce=0.230
#> n_total=36 (n1=12, n2=24): power=0.474, type1=0.027, ce=0.224
#> n_total=37 (n1=12, n2=25): power=0.479, type1=0.028, ce=0.255
#> n_total=38 (n1=13, n2=25): power=0.495, type1=0.030, ce=0.263
#> n_total=39 (n1=13, n2=26): power=0.492, type1=0.028, ce=0.257
#> n_total=40 (n1=13, n2=27): power=0.495, type1=0.028, ce=0.252
#> n_total=41 (n1=14, n2=27): power=0.500, type1=0.027, ce=0.288
#> n_total=42 (n1=14, n2=28): power=0.506, type1=0.028, ce=0.282
#> n_total=43 (n1=14, n2=29): power=0.511, type1=0.029, ce=0.316
#> n_total=44 (n1=15, n2=29): power=0.504, type1=0.024, ce=0.306
#> n_total=45 (n1=15, n2=30): power=0.516, type1=0.027, ce=0.300
#> n_total=46 (n1=15, n2=31): power=0.516, type1=0.026, ce=0.307
#> n_total=47 (n1=16, n2=31): power=0.522, type1=0.025, ce=0.321
#> n_total=48 (n1=16, n2=32): power=0.517, type1=0.023, ce=0.336
#> n_total=49 (n1=16, n2=33): power=0.522, type1=0.024, ce=0.341
#> n_total=50 (n1=17, n2=33): power=0.533, type1=0.025, ce=0.361
#> n_total=51 (n1=17, n2=34): power=0.530, type1=0.024, ce=0.355
#> n_total=52 (n1=17, n2=35): power=0.528, type1=0.022, ce=0.368
#> n_total=53 (n1=18, n2=35): power=0.535, type1=0.022, ce=0.376
#> n_total=54 (n1=18, n2=36): power=0.543, type1=0.024, ce=0.378
#> n_total=55 (n1=18, n2=37): power=0.543, type1=0.024, ce=0.373
#> n_total=56 (n1=19, n2=37): power=0.547, type1=0.023, ce=0.419
#> n_total=57 (n1=19, n2=38): power=0.546, type1=0.022, ce=0.428
#> n_total=58 (n1=19, n2=39): power=0.557, type1=0.025, ce=0.422
#> n_total=59 (n1=20, n2=39): power=0.557, type1=0.023, ce=0.470
#> n_total=60 (n1=20, n2=40): power=0.557, type1=0.022, ce=0.473
#> n_total=61 (n1=20, n2=41): power=0.556, type1=0.021, ce=0.469
#> n_total=62 (n1=21, n2=41): power=0.563, type1=0.022, ce=0.492
#> n_total=63 (n1=21, n2=42): power=0.573, type1=0.025, ce=0.497
#> n_total=64 (n1=21, n2=43): power=0.564, type1=0.021, ce=0.508
#> n_total=65 (n1=22, n2=43): power=0.573, type1=0.022, ce=0.516
#> n_total=66 (n1=22, n2=44): power=0.572, type1=0.022, ce=0.523
#> n_total=67 (n1=22, n2=45): power=0.578, type1=0.023, ce=0.536
#> n_total=68 (n1=23, n2=45): power=0.580, type1=0.022, ce=0.535
#> n_total=69 (n1=23, n2=46): power=0.580, type1=0.021, ce=0.550
#> n_total=70 (n1=23, n2=47): power=0.578, type1=0.021, ce=0.554
#> n_total=71 (n1=24, n2=47): power=0.587, type1=0.022, ce=0.567
#> n_total=72 (n1=24, n2=48): power=0.591, type1=0.023, ce=0.563
#> n_total=73 (n1=24, n2=49): power=0.586, type1=0.021, ce=0.573
#> n_total=74 (n1=25, n2=49): power=0.592, type1=0.021, ce=0.585
#> n_total=75 (n1=25, n2=50): power=0.591, type1=0.020, ce=0.592
#> n_total=76 (n1=25, n2=51): power=0.596, type1=0.022, ce=0.576
#> n_total=77 (n1=26, n2=51): power=0.597, type1=0.020, ce=0.596
#> n_total=78 (n1=26, n2=52): power=0.597, type1=0.020, ce=0.600
#> n_total=79 (n1=26, n2=53): power=0.601, type1=0.021, ce=0.591
#> n_total=80 (n1=27, n2=53): power=0.602, type1=0.020, ce=0.613
#> n_total=81 (n1=27, n2=54): power=0.608, type1=0.022, ce=0.606
#> n_total=82 (n1=27, n2=55): power=0.605, type1=0.020, ce=0.622
#> n_total=83 (n1=28, n2=55): power=0.611, type1=0.021, ce=0.620
#> n_total=84 (n1=28, n2=56): power=0.609, type1=0.020, ce=0.615
#> n_total=85 (n1=28, n2=57): power=0.611, type1=0.020, ce=0.632
#> n_total=86 (n1=29, n2=57): power=0.617, type1=0.021, ce=0.618
#> n_total=87 (n1=29, n2=58): power=0.614, type1=0.019, ce=0.643
#> n_total=88 (n1=29, n2=59): power=0.613, type1=0.019, ce=0.633
#> n_total=89 (n1=30, n2=59): power=0.619, type1=0.020, ce=0.626
#> n_total=90 (n1=30, n2=60): power=0.622, type1=0.020, ce=0.642
#> n_total=91 (n1=30, n2=61): power=0.620, type1=0.019, ce=0.644
#> n_total=92 (n1=31, n2=61): power=0.623, type1=0.019, ce=0.655
#> n_total=93 (n1=31, n2=62): power=0.623, type1=0.019, ce=0.643
#> n_total=94 (n1=31, n2=63): power=0.627, type1=0.020, ce=0.655
#> n_total=95 (n1=32, n2=63): power=0.631, type1=0.020, ce=0.669
#> n_total=96 (n1=32, n2=64): power=0.628, type1=0.019, ce=0.658
#> n_total=97 (n1=32, n2=65): power=0.632, type1=0.020, ce=0.654
#> n_total=98 (n1=33, n2=65): power=0.632, type1=0.019, ce=0.667
#> n_total=99 (n1=33, n2=66): power=0.634, type1=0.019, ce=0.675
#> n_total=100 (n1=33, n2=67): power=0.634, type1=0.019, ce=0.658
#> n_total=101 (n1=34, n2=67): power=0.635, type1=0.018, ce=0.668
#> n_total=102 (n1=34, n2=68): power=0.641, type1=0.020, ce=0.674
#> n_total=103 (n1=34, n2=69): power=0.638, type1=0.019, ce=0.677
#> n_total=104 (n1=35, n2=69): power=0.640, type1=0.018, ce=0.668
#> n_total=105 (n1=35, n2=70): power=0.641, type1=0.018, ce=0.682
#> n_total=106 (n1=35, n2=71): power=0.643, type1=0.019, ce=0.693
#> n_total=107 (n1=36, n2=71): power=0.645, type1=0.018, ce=0.698
#> n_total=108 (n1=36, n2=72): power=0.644, type1=0.018, ce=0.685
#> n_total=109 (n1=36, n2=73): power=0.646, type1=0.018, ce=0.699
#> n_total=110 (n1=37, n2=73): power=0.647, type1=0.017, ce=0.701
#> n_total=111 (n1=37, n2=74): power=0.648, type1=0.018, ce=0.685
#> n_total=112 (n1=37, n2=75): power=0.648, type1=0.017, ce=0.693
#> n_total=113 (n1=38, n2=75): power=0.650, type1=0.017, ce=0.705
#> n_total=114 (n1=38, n2=76): power=0.651, type1=0.017, ce=0.715
#> n_total=115 (n1=38, n2=77): power=0.652, type1=0.017, ce=0.688
#> n_total=116 (n1=39, n2=77): power=0.655, type1=0.018, ce=0.703
#> n_total=117 (n1=39, n2=78): power=0.653, type1=0.017, ce=0.712
#> n_total=118 (n1=39, n2=79): power=0.655, type1=0.017, ce=0.698
#> n_total=119 (n1=40, n2=79): power=0.660, type1=0.018, ce=0.709
#> n_total=120 (n1=40, n2=80): power=0.658, type1=0.017, ce=0.711
#> n_total=121 (n1=40, n2=81): power=0.657, type1=0.017, ce=0.721
#> n_total=122 (n1=41, n2=81): power=0.658, type1=0.016, ce=0.714
#> n_total=123 (n1=41, n2=82): power=0.662, type1=0.017, ce=0.709
#> n_total=124 (n1=41, n2=83): power=0.662, type1=0.017, ce=0.732
#> n_total=125 (n1=42, n2=83): power=0.663, type1=0.016, ce=0.723
#> n_total=126 (n1=42, n2=84): power=0.664, type1=0.016, ce=0.719
#> n_total=127 (n1=42, n2=85): power=0.665, type1=0.017, ce=0.729
#> n_total=128 (n1=43, n2=85): power=0.666, type1=0.016, ce=0.729
#> n_total=129 (n1=43, n2=86): power=0.666, type1=0.016, ce=0.723
#> n_total=130 (n1=43, n2=87): power=0.669, type1=0.017, ce=0.724
#> n_total=131 (n1=44, n2=87): power=0.671, type1=0.017, ce=0.734
#> n_total=132 (n1=44, n2=88): power=0.670, type1=0.016, ce=0.742
#> n_total=133 (n1=44, n2=89): power=0.671, type1=0.016, ce=0.720
#> n_total=134 (n1=45, n2=89): power=0.672, type1=0.016, ce=0.738
#> n_total=135 (n1=45, n2=90): power=0.673, type1=0.016, ce=0.741
#> n_total=136 (n1=45, n2=91): power=0.672, type1=0.016, ce=0.728
#> n_total=137 (n1=46, n2=91): power=0.675, type1=0.016, ce=0.740
#> n_total=138 (n1=46, n2=92): power=0.677, type1=0.016, ce=0.738
#> n_total=139 (n1=46, n2=93): power=0.675, type1=0.016, ce=0.754
#> n_total=140 (n1=47, n2=93): power=0.678, type1=0.016, ce=0.744
#> n_total=141 (n1=47, n2=94): power=0.679, type1=0.016, ce=0.737
#> n_total=142 (n1=47, n2=95): power=0.679, type1=0.016, ce=0.755
#> n_total=143 (n1=48, n2=95): power=0.682, type1=0.016, ce=0.748
#> n_total=144 (n1=48, n2=96): power=0.681, type1=0.016, ce=0.744
#> n_total=145 (n1=48, n2=97): power=0.681, type1=0.016, ce=0.753
#> n_total=146 (n1=49, n2=97): power=0.683, type1=0.015, ce=0.748
#> n_total=147 (n1=49, n2=98): power=0.683, type1=0.015, ce=0.758
#> n_total=148 (n1=49, n2=99): power=0.683, type1=0.015, ce=0.749
#> n_total=149 (n1=50, n2=99): power=0.686, type1=0.015, ce=0.755
#> n_total=150 (n1=50, n2=100): power=0.687, type1=0.016, ce=0.764
#> 
#> One-stage two-arm Bayes factor design
#> ------------------------------------
#> Mode: optimal
#> Status: No feasible one-stage two-arm design found.
#> Calibration: Bayesian
#> Optional freq. Type-I reporting: off
#> Design: no feasible design found
# }

# BF+0 directional test with plot
# \donttest{
ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.9,
               test = "BF+0", output = "plot")
#> n_total=10 (n1=5, n2=5): power=0.327, type1=0.029, ce=0.385
#> n_total=11 (n1=6, n2=5): power=0.372, type1=0.040, ce=0.343
#> n_total=12 (n1=6, n2=6): power=0.399, type1=0.040, ce=0.362
#> n_total=13 (n1=6, n2=7): power=0.386, type1=0.033, ce=0.367
#> n_total=14 (n1=7, n2=7): power=0.458, type1=0.052, ce=0.401
#> n_total=15 (n1=8, n2=7): power=0.437, type1=0.040, ce=0.363
#> n_total=16 (n1=8, n2=8): power=0.460, type1=0.046, ce=0.385
#> n_total=17 (n1=8, n2=9): power=0.480, type1=0.048, ce=0.387
#> n_total=18 (n1=9, n2=9): power=0.454, type1=0.036, ce=0.374
#> n_total=19 (n1=10, n2=9): power=0.500, type1=0.049, ce=0.440
#> n_total=20 (n1=10, n2=10): power=0.489, type1=0.041, ce=0.449
#> n_total=21 (n1=10, n2=11): power=0.479, type1=0.035, ce=0.433
#> n_total=22 (n1=11, n2=11): power=0.494, type1=0.037, ce=0.436
#> n_total=23 (n1=12, n2=11): power=0.495, type1=0.034, ce=0.464
#> n_total=24 (n1=12, n2=12): power=0.526, type1=0.042, ce=0.450
#> n_total=25 (n1=12, n2=13): power=0.511, type1=0.035, ce=0.472
#> n_total=26 (n1=13, n2=13): power=0.516, type1=0.034, ce=0.461
#> n_total=27 (n1=14, n2=13): power=0.528, type1=0.036, ce=0.477
#> n_total=28 (n1=14, n2=14): power=0.520, type1=0.031, ce=0.498
#> n_total=29 (n1=14, n2=15): power=0.545, type1=0.037, ce=0.505
#> n_total=30 (n1=15, n2=15): power=0.542, type1=0.034, ce=0.512
#> n_total=31 (n1=16, n2=15): power=0.533, type1=0.029, ce=0.497
#> n_total=32 (n1=16, n2=16): power=0.549, type1=0.032, ce=0.529
#> n_total=33 (n1=16, n2=17): power=0.552, type1=0.032, ce=0.523
#> n_total=34 (n1=17, n2=17): power=0.570, type1=0.036, ce=0.552
#> n_total=35 (n1=18, n2=17): power=0.558, type1=0.030, ce=0.540
#> n_total=36 (n1=18, n2=18): power=0.567, type1=0.032, ce=0.563
#> n_total=37 (n1=18, n2=19): power=0.575, type1=0.032, ce=0.551
#> n_total=38 (n1=19, n2=19): power=0.572, type1=0.030, ce=0.557
#> n_total=39 (n1=20, n2=19): power=0.587, type1=0.033, ce=0.580
#> n_total=40 (n1=20, n2=20): power=0.577, type1=0.029, ce=0.565
#> n_total=41 (n1=20, n2=21): power=0.581, type1=0.029, ce=0.594
#> n_total=42 (n1=21, n2=21): power=0.592, type1=0.031, ce=0.578
#> n_total=43 (n1=22, n2=21): power=0.586, type1=0.028, ce=0.604
#> n_total=44 (n1=22, n2=22): power=0.598, type1=0.030, ce=0.584
#> n_total=45 (n1=22, n2=23): power=0.591, type1=0.027, ce=0.607
#> n_total=46 (n1=23, n2=23): power=0.604, type1=0.030, ce=0.602
#> n_total=47 (n1=24, n2=23): power=0.604, type1=0.029, ce=0.608
#> n_total=48 (n1=24, n2=24): power=0.602, type1=0.027, ce=0.606
#> n_total=49 (n1=24, n2=25): power=0.607, type1=0.027, ce=0.616
#> n_total=50 (n1=25, n2=25): power=0.607, type1=0.026, ce=0.617
#> n_total=51 (n1=26, n2=25): power=0.615, type1=0.028, ce=0.621
#> n_total=52 (n1=26, n2=26): power=0.612, type1=0.026, ce=0.642
#> n_total=53 (n1=26, n2=27): power=0.619, type1=0.027, ce=0.621
#> n_total=54 (n1=27, n2=27): power=0.618, type1=0.026, ce=0.639
#> n_total=55 (n1=28, n2=27): power=0.616, type1=0.024, ce=0.630
#> n_total=56 (n1=28, n2=28): power=0.623, type1=0.025, ce=0.644
#> n_total=57 (n1=28, n2=29): power=0.621, type1=0.024, ce=0.633
#> n_total=58 (n1=29, n2=29): power=0.633, type1=0.027, ce=0.649
#> n_total=59 (n1=30, n2=29): power=0.630, type1=0.025, ce=0.640
#> n_total=60 (n1=30, n2=30): power=0.631, type1=0.024, ce=0.646
#> n_total=61 (n1=30, n2=31): power=0.633, type1=0.024, ce=0.642
#> n_total=62 (n1=31, n2=31): power=0.631, type1=0.023, ce=0.650
#> n_total=63 (n1=32, n2=31): power=0.642, type1=0.025, ce=0.654
#> n_total=64 (n1=32, n2=32): power=0.635, type1=0.022, ce=0.653
#> n_total=65 (n1=32, n2=33): power=0.647, type1=0.025, ce=0.664
#> n_total=66 (n1=33, n2=33): power=0.643, type1=0.023, ce=0.656
#> n_total=67 (n1=34, n2=33): power=0.645, type1=0.023, ce=0.670
#> n_total=68 (n1=34, n2=34): power=0.648, type1=0.023, ce=0.658
#> n_total=69 (n1=34, n2=35): power=0.647, type1=0.023, ce=0.670
#> n_total=70 (n1=35, n2=35): power=0.655, type1=0.024, ce=0.660
#> n_total=71 (n1=36, n2=35): power=0.651, type1=0.022, ce=0.671
#> n_total=72 (n1=36, n2=36): power=0.657, type1=0.024, ce=0.662
#> n_total=73 (n1=36, n2=37): power=0.655, type1=0.022, ce=0.675
#> n_total=74 (n1=37, n2=37): power=0.662, type1=0.024, ce=0.668
#> n_total=75 (n1=38, n2=37): power=0.660, type1=0.022, ce=0.678
#> n_total=76 (n1=38, n2=38): power=0.660, type1=0.022, ce=0.681
#> n_total=77 (n1=38, n2=39): power=0.665, type1=0.023, ce=0.680
#> n_total=78 (n1=39, n2=39): power=0.663, type1=0.022, ce=0.698
#> n_total=79 (n1=40, n2=39): power=0.670, type1=0.023, ce=0.685
#> n_total=80 (n1=40, n2=40): power=0.665, type1=0.021, ce=0.701
#> n_total=81 (n1=40, n2=41): power=0.675, type1=0.024, ce=0.693
#> n_total=82 (n1=41, n2=41): power=0.671, type1=0.022, ce=0.704
#> n_total=83 (n1=42, n2=41): power=0.672, type1=0.021, ce=0.699
#> n_total=84 (n1=42, n2=42): power=0.674, type1=0.022, ce=0.712
#> n_total=85 (n1=42, n2=43): power=0.674, type1=0.021, ce=0.707
#> n_total=86 (n1=43, n2=43): power=0.680, type1=0.022, ce=0.714
#> n_total=87 (n1=44, n2=43): power=0.674, type1=0.020, ce=0.712
#> n_total=88 (n1=44, n2=44): power=0.684, type1=0.023, ce=0.716
#> n_total=89 (n1=44, n2=45): power=0.679, type1=0.020, ce=0.714
#> n_total=90 (n1=45, n2=45): power=0.685, type1=0.022, ce=0.717
#> n_total=91 (n1=46, n2=45): power=0.682, type1=0.020, ce=0.722
#> n_total=92 (n1=46, n2=46): power=0.682, type1=0.020, ce=0.719
#> n_total=93 (n1=46, n2=47): power=0.687, type1=0.021, ce=0.735
#> n_total=94 (n1=47, n2=47): power=0.684, type1=0.019, ce=0.720
#> n_total=95 (n1=48, n2=47): power=0.690, type1=0.021, ce=0.735
#> n_total=96 (n1=48, n2=48): power=0.685, type1=0.019, ce=0.724
#> n_total=97 (n1=48, n2=49): power=0.693, type1=0.021, ce=0.734
#> n_total=98 (n1=49, n2=49): power=0.688, type1=0.019, ce=0.724
#> n_total=99 (n1=50, n2=49): power=0.691, type1=0.020, ce=0.734
#> n_total=100 (n1=50, n2=50): power=0.693, type1=0.020, ce=0.728
#> n_total=101 (n1=50, n2=51): power=0.692, type1=0.019, ce=0.733
#> n_total=102 (n1=51, n2=51): power=0.696, type1=0.020, ce=0.731
#> n_total=103 (n1=52, n2=51): power=0.693, type1=0.018, ce=0.733
#> n_total=104 (n1=52, n2=52): power=0.699, type1=0.020, ce=0.734
#> n_total=105 (n1=52, n2=53): power=0.696, type1=0.018, ce=0.732
#> n_total=106 (n1=53, n2=53): power=0.701, type1=0.020, ce=0.745
#> n_total=107 (n1=54, n2=53): power=0.700, type1=0.019, ce=0.731
#> n_total=108 (n1=54, n2=54): power=0.701, type1=0.019, ce=0.746
#> n_total=109 (n1=54, n2=55): power=0.703, type1=0.019, ce=0.732
#> n_total=110 (n1=55, n2=55): power=0.702, type1=0.018, ce=0.744
#> n_total=111 (n1=56, n2=55): power=0.704, type1=0.019, ce=0.731
#> n_total=112 (n1=56, n2=56): power=0.704, type1=0.018, ce=0.742
#> n_total=113 (n1=56, n2=57): power=0.708, type1=0.019, ce=0.735
#> n_total=114 (n1=57, n2=57): power=0.706, type1=0.018, ce=0.740
#> n_total=115 (n1=58, n2=57): power=0.708, type1=0.018, ce=0.739
#> n_total=116 (n1=58, n2=58): power=0.709, type1=0.018, ce=0.739
#> n_total=117 (n1=58, n2=59): power=0.709, type1=0.018, ce=0.746
#> n_total=118 (n1=59, n2=59): power=0.711, type1=0.018, ce=0.740
#> n_total=119 (n1=60, n2=59): power=0.710, type1=0.017, ce=0.748
#> n_total=120 (n1=60, n2=60): power=0.713, type1=0.018, ce=0.751
#> n_total=121 (n1=60, n2=61): power=0.713, type1=0.018, ce=0.755
#> n_total=122 (n1=61, n2=61): power=0.716, type1=0.018, ce=0.754
#> n_total=123 (n1=62, n2=61): power=0.714, type1=0.017, ce=0.768
#> n_total=124 (n1=62, n2=62): power=0.718, type1=0.018, ce=0.757
#> n_total=125 (n1=62, n2=63): power=0.716, type1=0.017, ce=0.770
#> n_total=126 (n1=63, n2=63): power=0.717, type1=0.017, ce=0.760
#> n_total=127 (n1=64, n2=63): power=0.719, type1=0.017, ce=0.770
#> n_total=128 (n1=64, n2=64): power=0.718, type1=0.017, ce=0.762
#> n_total=129 (n1=64, n2=65): power=0.721, type1=0.017, ce=0.771
#> n_total=130 (n1=65, n2=65): power=0.720, type1=0.017, ce=0.766
#> n_total=131 (n1=66, n2=65): power=0.724, type1=0.018, ce=0.771
#> n_total=132 (n1=66, n2=66): power=0.721, type1=0.016, ce=0.768
#> n_total=133 (n1=66, n2=67): power=0.725, type1=0.017, ce=0.771
#> n_total=134 (n1=67, n2=67): power=0.724, type1=0.017, ce=0.771
#> n_total=135 (n1=68, n2=67): power=0.724, type1=0.017, ce=0.772
#> n_total=136 (n1=68, n2=68): power=0.726, type1=0.017, ce=0.776
#> n_total=137 (n1=68, n2=69): power=0.726, type1=0.017, ce=0.771
#> n_total=138 (n1=69, n2=69): power=0.728, type1=0.017, ce=0.786
#> n_total=139 (n1=70, n2=69): power=0.727, type1=0.016, ce=0.772
#> n_total=140 (n1=70, n2=70): power=0.731, type1=0.017, ce=0.784
#> n_total=141 (n1=70, n2=71): power=0.728, type1=0.016, ce=0.773
#> n_total=142 (n1=71, n2=71): power=0.733, type1=0.017, ce=0.783
#> n_total=143 (n1=72, n2=71): power=0.730, type1=0.016, ce=0.773
#> n_total=144 (n1=72, n2=72): power=0.732, type1=0.016, ce=0.782
#> n_total=145 (n1=72, n2=73): power=0.732, type1=0.016, ce=0.773
#> n_total=146 (n1=73, n2=73): power=0.733, type1=0.016, ce=0.780
#> n_total=147 (n1=74, n2=73): power=0.735, type1=0.017, ce=0.773
#> n_total=148 (n1=74, n2=74): power=0.733, type1=0.016, ce=0.779
#> n_total=149 (n1=74, n2=75): power=0.736, type1=0.016, ce=0.773
#> n_total=150 (n1=75, n2=75): power=0.734, type1=0.015, ce=0.777

# }
```
