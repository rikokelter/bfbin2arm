# Bayes factor BF+0 for the directional alternative vs point-null

Computes the Bayes factor \\BF\_{+0}\\ comparing the directional
alternative hypothesis \\H\_+\\ (p_2 \> p_1) against the point-null
\\H_0\\ (p_1 = p_2).

## Usage

``` r
twoarmbinbf_plus0_direct(
  y1,
  y2,
  n1,
  n2,
  a_0_a = 1,
  b_0_a = 1,
  a_1_a = 1,
  b_1_a = 1,
  a_2_a = 1,
  b_2_a = 1
)
```

## Arguments

- y1, y2:

  Integer counts of successes in arms 1 and 2.

- n1, n2:

  Integer sample sizes in arms 1 and 2.

- a_0_a, b_0_a:

  Shape parameters of the **analysis** prior for the common response
  probability under \\H_0\\.

- a_1_a, b_1_a:

  Shape parameters of the **analysis** prior for the response
  probability in arm 1 under \\H\_+\\.

- a_2_a, b_2_a:

  Shape parameters of the **analysis** prior for the response
  probability in arm 2 under \\H\_+\\.

## Value

Numeric scalar; the Bayes factor \\BF\_{+0} = m\_+(y_1, y_2) / m_0(y_1,
y_2)\\.

## Details

Both marginal likelihoods are formed using the **analysis** priors,
which represent inferential beliefs at the time the data are evaluated.
The design priors are used only for computing Bayesian operating
characteristics (predictive power / type-I error) and play no role here.
