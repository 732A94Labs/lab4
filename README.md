# lab4 + Bonus lab

<!-- badges: start -->
[![R-CMD-check](https://github.com/732A94Labs/lab4/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/732A94Labs/lab4/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`lab4` is the Advanced R Programming Lab 4 package. It delivers a compact linear regression toolkit: ordinary least squares core computations, optional QR-based solver, model methods for summaries/predictions/plots, and a Linköping University `ggplot2` theme.

## Install

```r
# install.packages("pak")
pak::pak("732A94Labs/lab4")
```

## Quick start

```r
library(lab4)

fit <- linreg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)
fit$summary()  # coefficients, t-stats, sigma
pl <- fit$plot()     # residual diagnostics
pl  # to display the generated plot
```

## Linköping theme

```r
library(ggplot2)

ggplot(iris, aes(Sepal.Length, Petal.Length, colour = Species)) +
  geom_point() +
  liu_theme()
```

## What's inside

- `linreg()` fits with coefficients, residuals, predictions, and summaries.
- Unit tests and a vignette covering a minimal regression workflow.
- Linköping University styling via a custom `ggplot2` theme.
