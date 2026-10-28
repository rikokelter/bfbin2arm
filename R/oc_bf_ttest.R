# ============================================================
# Decision boundaries and operating characteristics
# ============================================================

p_compelling_h1_one_sample <- function(delta, n, c1, c2) {
  df <- n - 1
  ncp <- sqrt(n) * delta
  
  stats::pt(c1, df = df, ncp = ncp) +
    stats::pt(
      c2,
      df = df,
      ncp = ncp,
      lower.tail = FALSE
    )
}


p_compelling_h1_two_sample_pooled <- function(
    delta,
    n1,
    n2,
    c1,
    c2
) {
  df <- n1 + n2 - 2
  n_eff <- n_eff_two_sample(n1, n2)
  ncp <- sqrt(n_eff) * delta
  
  stats::pt(c1, df = df, ncp = ncp) +
    stats::pt(
      c2,
      df = df,
      ncp = ncp,
      lower.tail = FALSE
    )
}


p_compelling_h0_interval <- function(
    delta,
    df,
    n_eff,
    a0,
    b0
) {
  ncp <- sqrt(n_eff) * delta
  
  stats::pt(b0, df = df, ncp = ncp) -
    stats::pt(a0, df = df, ncp = ncp)
}

# Same for a pooled two-sample t-test.
p_compelling_h1_two_sample_pooled_upper <- function(delta, n1, n2, c1) {
  df <- n1 + n2 - 2L
  n_eff <- n_eff_two_sample(n1, n2)
  ncp <- sqrt(n_eff) * delta
  pnt_stable(c1, df = df, ncp = ncp, lower.tail = FALSE)
}

# Generic probability of a union of t-statistic intervals.
# Each row of `intervals` is c(lower, upper).
p_t_in_intervals <- function(delta, df, n_eff, intervals) {
  if (is.null(intervals) || !length(intervals)) return(0)
  
  intervals <- as.matrix(intervals)
  if (ncol(intervals) != 2L) {
    stop("`intervals` must have two columns: lower and upper.", call. = FALSE)
  }
  
  cdf_nct <- function(x, df, ncp) {
    if (is.infinite(x) && x < 0) return(0)
    if (is.infinite(x) && x > 0) return(1)
    pnt_stable(x, df = df, ncp = ncp)
  }
  
  ncp <- sqrt(n_eff) * delta
  
  sum(vapply(
    seq_len(nrow(intervals)),
    function(i) {
      lower <- intervals[i, 1L]
      upper <- intervals[i, 2L]
      cdf_nct(upper, df = df, ncp = ncp) -
        cdf_nct(lower, df = df, ncp = ncp)
    },
    numeric(1)
  ))
}
