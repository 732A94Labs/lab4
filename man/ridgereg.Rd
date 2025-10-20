% Generated manually because roxygen2 is unavailable in the execution environment
\name{ridgereg}
\alias{ridgereg}
\alias{print.ridgereg}
\alias{coef.ridgereg}
\alias{residuals.ridgereg}
\alias{fitted.ridgereg}
\alias{predict.ridgereg}
\alias{summary.ridgereg}
\alias{print.summary.ridgereg}
\title{Ridge Regression via Closed-Form Solution}
\usage{
ridgereg(formula, data, lambda = 0)
}
\arguments{
  \item{formula}{A formula specifying the response and predictors.}
  \item{data}{A data frame (or tibble) containing the variables in the model.}
  \item{lambda}{A non-negative numeric value controlling the amount of shrinkage.}
}
\value{
An object of class \code{ridgereg} containing fitted coefficients, residuals
and summary information for the ridge regression model.
}
\description{
Fits a ridge regression model using the closed-form normal equations applied to
centre- and scale-standardised covariates. Convenience methods are provided for
printing, extracting coefficients, computing fitted values, predicting on new
observations and summarising the model fit.
}
\details{
The function normalises all covariates to mean zero and unit variance prior to
fitting, and then back-transforms the resulting estimates so that the reported
coefficients correspond to the original data scale. The intercept is not
penalised. Predictions on new data rely on the stored design information from
the training set.
}
\examples{
if (requireNamespace("datasets", quietly = TRUE)) {
  data(iris)
  mod <- ridgereg(Petal.Length ~ Sepal.Width + Sepal.Length, data = iris, lambda = 1)
  coef(mod)
  head(predict(mod))
}
}
