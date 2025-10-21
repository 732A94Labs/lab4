#' Ridge Regression via Closed-Form Solution
#'
#' Fits a ridge regression model using normal-equation algebra with centred and
#' scaled covariates. The function stores the fitted coefficients on the original
#' data scale together with the key model diagnostics and supports standard S3
#' methods for printing, prediction and summarising the result.
#'
#' @param formula A model formula specifying the response and predictors.
#' @param data A `data.frame` (or tibble) containing the model variables.
#' @param lambda A non-negative numeric penalty parameter controlling the shrinkage.
#'
#' @return An object of class `ridgereg`.
#' @export
#'
#' @importFrom stats model.frame model.matrix model.response delete.response pt
ridgereg <- function(formula, data, lambda = 0) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("`formula` must be provided and must be a formula.", call. = FALSE)
  }
  if (missing(data)) {
    stop("`data` must be supplied.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("`lambda` must be a single non-negative numeric value.", call. = FALSE)
  }

  cl <- match.call()

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  if (!is.numeric(y)) {
    stop("The response must be numeric for ridge regression.", call. = FALSE)
  }

  terms <- attr(mf, "terms")
  X <- stats::model.matrix(terms, mf)
  contrasts <- attr(X, "contrasts")
  xlevels <- stats::.getXlevels(terms, mf)

  n <- nrow(X)
  p <- ncol(X)

  if (n == 0 || p == 0) {
    stop("Model matrix is empty; please check the supplied data.", call. = FALSE)
  }

  # Extract covariate matrix without the intercept for scaling purposes.
  has_intercept <- colnames(X)[1] == "(Intercept)"
  if (has_intercept) {
    X_cov <- X[, -1, drop = FALSE]
  } else {
    X_cov <- X
  }

  center <- if (ncol(X_cov) > 0) {
    colMeans(X_cov)
  } else {
    numeric(0)
  }
  scale <- if (ncol(X_cov) > 0) {
    out <- apply(X_cov, 2, stats::sd)
    out[out == 0 | is.na(out)] <- 1
    out
  } else {
    numeric(0)
  }

  if (ncol(X_cov) > 0) {
    X_scaled_cov <- sweep(sweep(X_cov, 2, center, FUN = "-"), 2, scale, FUN = "/")
  } else {
    X_scaled_cov <- X_cov
  }

  if (has_intercept) {
    X_scaled <- cbind(`(Intercept)` = 1, X_scaled_cov)
  } else {
    X_scaled <- X_scaled_cov
  }

  penalty <- diag(ncol(X_scaled))
  if (has_intercept) {
    penalty[1, 1] <- 0
  }

  XtX <- crossprod(X_scaled)
  Xty <- crossprod(X_scaled, y)
  ridge_inv <- tryCatch(
    solve(XtX + lambda * penalty),
    error = function(e) {
      stop("Failed to invert penalised cross-product matrix; check for collinearity.", call. = FALSE)
    }
  )
  sol <- ridge_inv %*% Xty

  coef_scaled <- as.numeric(sol)
  names(coef_scaled) <- colnames(X_scaled)

  if (has_intercept) {
    coef_original <- numeric(length = ncol(X))
    names(coef_original) <- colnames(X)

    coef_original[1] <- coef_scaled[1]
    if (length(coef_scaled) > 1) {
      covar_names <- colnames(X)[-1]
      adj <- coef_scaled[-1] / scale[covar_names]
      coef_original[-1] <- adj
      coef_original[1] <- coef_original[1] - sum(adj * center[covar_names])
    }
  } else {
    coef_original <- coef_scaled / scale[colnames(X)]
    names(coef_original) <- colnames(X)
  }

  fitted <- as.numeric(X %*% coef_original)
  residuals <- y - fitted

  hat_diag <- rowSums((X_scaled %*% ridge_inv) * X_scaled)
  effective_df <- sum(hat_diag)
  df_residual <- max(n - effective_df, .Machine$double.eps)

  rss <- sum(residuals^2)
  sigma2 <- rss / df_residual
  sigma <- sqrt(sigma2)

  vcov_scaled <- sigma2 * ridge_inv %*% XtX %*% ridge_inv

  trans <- matrix(0, nrow = length(coef_original), ncol = length(coef_scaled))
  rownames(trans) <- names(coef_original)
  colnames(trans) <- names(coef_scaled)

  if (has_intercept) {
    trans[1, 1] <- 1
    if (length(coef_scaled) > 1) {
      covar_names <- colnames(X_scaled)[-1]
      trans[1, covar_names] <- -center[covar_names] / scale[covar_names]
      for (nm in covar_names) {
        trans[nm, nm] <- 1 / scale[nm]
      }
    }
  } else if (length(coef_scaled) > 0) {
    covar_names <- colnames(X_scaled)
    for (nm in covar_names) {
      trans[nm, nm] <- 1 / scale[nm]
    }
  }

  vcov_original <- trans %*% vcov_scaled %*% t(trans)

  se <- sqrt(diag(vcov_original))
  t_values <- coef_original / se
  p_values <- 2 * stats::pt(-abs(t_values), df = df_residual)

  summary_table <- data.frame(
    Estimate = coef_original,
    `Std. Error` = se,
    `t value` = t_values,
    `Pr(>|t|)` = p_values,
    check.names = FALSE
  )

  out <- list(
    call = cl,
    formula = formula,
    coefficients = coef_original,
    coefficients_scaled = coef_scaled,
    fitted.values = fitted,
    residuals = residuals,
    lambda = lambda,
    center = center,
    scale = scale,
    sigma = sigma,
    df_residual = df_residual,
    effective_df = effective_df,
    vcov = vcov_original,
    summary_table = summary_table,
    terms = terms,
    contrasts = contrasts,
    xlevels = xlevels
  )
  class(out) <- "ridgereg"
  out
}

#' @export
print.ridgereg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
coef.ridgereg <- function(object, ...) {
  object$coefficients
}

#' @export
residuals.ridgereg <- function(object, ...) {
  object$residuals
}

#' @export
fitted.ridgereg <- function(object, ...) {
  object$fitted.values
}

#' @export
predict.ridgereg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  Terms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(Terms, newdata, na.action = stats::na.pass, xlev = object$xlevels)
  X_new <- stats::model.matrix(Terms, mf, contrasts.arg = object$contrasts)

  missing_cols <- setdiff(names(object$coefficients), colnames(X_new))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      X_new <- cbind(X_new, 0)
      colnames(X_new)[ncol(X_new)] <- col
    }
  }
  X_new <- X_new[, names(object$coefficients), drop = FALSE]
  as.numeric(X_new %*% object$coefficients)
}

#' @export
summary.ridgereg <- function(object, ...) {
  res <- list(
    call = object$call,
    coefficients = object$summary_table,
    sigma = object$sigma,
    df_residual = object$df_residual,
    lambda = object$lambda
  )
  class(res) <- "summary.ridgereg"
  res
}

#' @export
print.summary.ridgereg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  cat(
    sprintf(
      "\nResidual standard error: %.4f on %.1f effective degrees of freedom\n",
      x$sigma,
      x$df_residual
    )
  )
  cat(sprintf("Lambda: %.4f\n", x$lambda))
  invisible(x)
}
