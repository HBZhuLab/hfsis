<!-- README.md is generated from README.Rmd. Please edit that file. -->

# hfsis

hfsis implements variable screening for synchronous, equally spaced
high-frequency regression increments. It provides marginal spot independence
screening (M-SIS) and partial-residual spot iterative screening (PR-SISIS).
M-SIS selects a local window from the marginal correlation vector. PR-SISIS
updates the selected window as the active set grows and computes correlation
columns only when they are required.

## Installation

```r
install.packages(
  "hfsis_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```
or

```r
devtools::install_github("HBZhuLab/hfsis")
```

## Quick start: M-SIS and PR-SISIS


``` r
library(hfsis)
data(hf_sis_example)

dat <- hf_sis_example$data
x <- as.matrix(dat[hf_sis_example$x_names])
y <- dat[[hf_sis_example$response]]

ctrl <- hfsis_control(
  windows = hf_sis_example$windows,
  truncation = "bipower"
)
inputs <- prepare_hfsis(
  x = x,
  y = y,
  target_index = hf_sis_example$target_index,
  control = ctrl
)

fit_m <- m_sis(
  inputs,
  list_size = 10L
)

fit_pr <- pr_sisis(
  inputs,
  model_cap = 10L,
  block_size = 2L
)

selected_variables(fit_m)
```

```
##  [1] "x03" "x04" "x05" "x19" "x08" "x18" "x20" "x13" "x24" "x23"
```

``` r
selected_variables(fit_pr)
```

```
##  [1] "x03" "x04" "x05" "x20" "x22" "x12" "x17" "x10" "x07" "x15"
```

``` r
compare_screens(fit_m, fit_pr)
```

```
## HF-SIS screening comparison
##   selected variables: M-SIS 10 / PR-SISIS 10 
##   overlap: 4 
##   M-SIS / PR-SISIS first windows: 128 vs 128 (same) 
##   PR-SISIS final window: 128
```

The prepared object stores the common truncated local data and marginal window
summaries. M-SIS uses the marginal-only window selector. PR-SISIS reselects the
window at each update after incorporating the correlation columns associated
with the current active set. Complete PR-SISIS fits provide a full ranking.
Exact score ties use the original column order.

The main diagnostic plots are:


``` r
plot(fit_m, type = "scores")
plot(fit_m, type = "window")
plot(fit_pr, type = "path")
plot(compare_screens(fit_m, fit_pr))
```

## Repeated screening


``` r
data(hf_sis_panel_example)
panel_fits <- unlist(lapply(hf_sis_panel_example$sessions, function(session) {
  prepared <- prepare_hfsis(
    as.matrix(session[hf_sis_panel_example$x_names]),
    session[[hf_sis_panel_example$response]],
    hf_sis_panel_example$target_index,
    hfsis_control(hf_sis_panel_example$windows)
  )
  list(
    m_sis(prepared, list_size = 5),
    pr_sisis(prepared, model_cap = 5L, block_size = 2L)
  )
}), recursive = FALSE)

frequencies <- selection_frequency(
  panel_fits,
  metadata = hf_sis_panel_example$metadata
)
head(frequencies)
```

```
##   variable method selected_count fit_count frequency mean_rank median_rank
## 1      x01  M-SIS              0         5         0      22.0          21
## 2      x02  M-SIS              0         5         0      20.6          21
## 3      x03  M-SIS              5         5         1       1.0           1
## 4      x04  M-SIS              5         5         1       2.2           2
## 5      x05  M-SIS              5         5         1       3.6           3
## 6      x06  M-SIS              0         5         0      21.6          22
##     cluster active
## 1 cluster_1   TRUE
## 2 cluster_1   TRUE
## 3 cluster_1   TRUE
## 4 cluster_1   TRUE
## 5 cluster_1   TRUE
## 6 cluster_2  FALSE
```

``` r
cluster_share(panel_fits, hf_sis_panel_example$metadata)
```

```
##      method   cluster variable_count selection_count fit_count mean_frequency
## 1     M-SIS cluster_1              5              15         5           0.60
## 2     M-SIS cluster_2              5               1         5           0.04
## 3     M-SIS cluster_3              5               2         5           0.08
## 4     M-SIS cluster_4              5               4         5           0.16
## 5     M-SIS cluster_5              5               3         5           0.12
## 6  PR-SISIS cluster_1              5              15         5           0.60
## 7  PR-SISIS cluster_2              5               3         5           0.12
## 8  PR-SISIS cluster_3              5               1         5           0.04
## 9  PR-SISIS cluster_4              5               2         5           0.08
## 10 PR-SISIS cluster_5              5               4         5           0.16
##    selection_share
## 1             0.60
## 2             0.04
## 3             0.08
## 4             0.16
## 5             0.12
## 6             0.60
## 7             0.12
## 8             0.04
## 9             0.08
## 10            0.16
```

Selection frequency is descriptive empirical recurrence. It is not a support
probability, p-value, or significance measure.

## Input conventions

- Rows are consecutive equally spaced increments or returns, not price levels.
- Columns of x are synchronous candidate covariates; y is the response
  increment.
- target_index is the first increment row in every forward candidate window.
- delta is observation spacing on the normalized time scale and defaults to
  1 / nrow(x).
- Column order is preserved and determines exact-tie resolution.
- Rows are never silently dropped, imputed, or standardized outside the
  correlation normalization defined by the method.

## Scope

Version 0.1.0 assumes synchronous, equally spaced increments and uses forward
local windows. The package does not adjust for asynchronous observations,
market microstructure noise, or price staleness.

## Citation

Please cite Haibin Zhu, “Ultra-high-dimensional spot support screening in
continuous-time regression,” together with hfsis version 0.1.0. The installed
package citation is available through:

```r
citation("hfsis")
```
