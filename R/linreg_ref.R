#' @title Linear Regression Class
#' @description The Class for fitting and analyzing ordinary least squares linear regression models.
#' @details This class uses simple linear algebra to compute regression coefficients, fitted values,
#' residuals, and some key statistical metrics alongside several methods for inspecting the fitted model.
#'
#' @field formula An object of class \link{formula} representing the model.
#' @field data A \code{data.frame} containing the variables in the model.
#' @field data_name A character string with the name of the data frame.
#' @field coefficients A named numeric vector of the regression coefficients.
#' @field fitted_values A numeric vector of the model's fitted (predicted) values.
#' @field residuals A numeric vector of the model's residuals.
#' @field df_residual A number representing the residual degrees of freedom.
#' @field sigma A number representing the residual standard error.
#' @field var_hat_beta A matrix representing the variance-covariance matrix of the regression coefficients.
#' @field summary_table A \code{data.frame} containing the estimates, standard errors, t-values, and p-values for the coefficients.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(formula, data)}}{
#'     Initializes a new linreg object.
#'     \strong{Arguments:}
#'     \itemize{
#'       \item \code{formula}: A model \link{formula}.
#'       \item \code{data}: A \code{data.frame} containing the data.
#'     }
#'     \strong{Returns:} A new \code{linreg} object.
#'   }
#'   \item{\code{print()}}{Prints the model call and the estimated coefficients.}
#'   \item{\code{plot()}}{
#'     Generates two diagnostic plots using \code{ggplot2}: a "Residuals vs Fitted" plot and a
#'     "Scale-Location" plot.
#'     \strong{Returns:} A \code{patchwork} object containing the combined plots, which is also printed to the screen.
#'   }
#'   \item{\code{pred()}}{Returns the vector of fitted values.}
#'   \item{\code{resid()}}{Returns the vector of residuals.}
#'   \item{\code{coef()}}{Returns the named vector of regression coefficients.}
#'   \item{\code{summary()}}{Prints a detailed summary of the model, including coefficient statistics and residual standard error.}
#' }
#'
#' @examples
#' # Create a new linear regression model using the iris dataset
#' linreg_model <- linreg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)
#'
#' # Print the model (this is called automatically when the object name is typed)
#' linreg_model
#'
#' # Get the model coefficients
#' linreg_model$coef()
#'
#' # Get the fitted values
#' linreg_model$pred()
#'
#' # Generate and view diagnostic plots
#' linreg_model$plot()
#'
#' # Print a detailed summary of the model
#' linreg_model$summary()
#'
#' @import methods
#' @import ggplot2
#' @import patchwork
#' @importFrom stats model.frame model.matrix model.response pt
#' @export linreg
#' @exportClass linreg
linreg <- setRefClass(
  "linreg",
  fields = list(
    formula = "formula",
    data = "data.frame",
    data_name = "character",
    coefficients = "numeric",
    fitted_values = "numeric",
    residuals = "numeric",
    df_residual = "numeric",
    sigma = "numeric",
    var_hat_beta = "matrix",
    summary_table = "data.frame"
  ),
  methods = list(
    initialize = function(formula, data, ...) {
      if (missing(formula) || !inherits(formula, "formula")) {
        stop("`formula` must be provided as a formula.")
      }
      if (missing(data) || !is.data.frame(data)) {
        stop("`data` must be provided as a data.frame.")
      }

      call <- sys.call(-1)
      data_label <- ""
      arg <- tryCatch(substitute(data), error = function(e) NULL)
      if (!is.null(arg) && is.symbol(arg)) {
        data_label <- deparse(arg)
      }

      mf <- model.frame(formula, data)
      y  <- model.response(mf)
      X  <- model.matrix(attr(mf, "terms"), mf)

      # beta_hat = (X'X)^(-1) X'y
      coef_hat <- as.numeric(solve(crossprod(X)) %*% crossprod(X, y))
      names(coef_hat) <- colnames(X)

      # fitted and residuals
      fitted_vals <- as.numeric(X %*% coef_hat)
      resid_vec   <- as.numeric(y - fitted_vals)

      # df = n - p
      df_resid_val <- length(y) - ncol(X)
      if (df_resid_val <= 0) stop("Residual degrees of freedom must be positive.")

      # sigma^2 = RSS / df
      rss          <- crossprod(resid_vec)
      sigma2_hat   <- as.numeric(rss / df_resid_val)

      # Var(beta_hat) = sigma^2 (X'X)^(-1)
      V_beta <- sigma2_hat * solve(crossprod(X))

      # t-stats and p-values
      t_values <- coef_hat / sqrt(diag(V_beta))
      p_values <- 2 * pt(-abs(t_values), df_resid_val)

      coef_table <- data.frame(
        Estimate   = coef_hat,
        `Std. Error` = sqrt(diag(V_beta)),
        `t value`    = t_values,
        `Pr(>|t|)`   = p_values,
        check.names = FALSE
      )

      # ---- assign to fields (no name shadowing now) ----
      .self$formula        <- formula
      .self$data           <- data
      .self$data_name      <- data_label
      .self$coefficients   <- coef_hat
      .self$fitted_values  <- fitted_vals
      .self$residuals      <- resid_vec
      .self$df_residual    <- df_resid_val
      .self$sigma          <- sqrt(sigma2_hat)
      .self$var_hat_beta   <- V_beta
      .self$summary_table  <- coef_table

      callSuper(...)
    },

    print = function(...) {
      cat(sprintf(
        "linreg(formula = %s, data = %s)\n",
        deparse(.self$formula),
        .self$data_name
      ))

      cat("Coefficients:\n")
      cat(paste(sprintf("%15s", names(.self$coefficients)), collapse = ""), "\n", sep = "")
      cat(paste(sprintf("%15s", format(.self$coefficients, digits = 4, trim = TRUE)), collapse = ""), "\n", sep = "")

      invisible(.self)
    },

    pred = function() {
      .self$fitted_values
    },

    resid = function() {
      .self$residuals
    },

    coef = function() {
      .self$coefficients
    },

    plot = function(...) {
      mf <- stats::model.frame(.self$formula, .self$data)
      X  <- stats::model.matrix(attr(mf, "terms"), mf)

      XtX_inv <- chol2inv(chol(crossprod(X)))
      h <- rowSums((X %*% XtX_inv) * X)

      rs   <- .self$residuals / (.self$sigma * sqrt(pmax(1e-12, 1 - h)))
      sl_y <- sqrt(abs(rs))

      df <- data.frame(
        fitted = .self$fitted_values,
        resid  = .self$residuals,
        sl_y   = sl_y
      )

      sm1 <- stats::lowess(df$fitted, df$resid,  f = 2/3)
      sm2 <- stats::lowess(df$fitted, df$sl_y,  f = 2/3)
      sm1df <- data.frame(x = sm1$x, y = sm1$y)
      sm2df <- data.frame(x = sm2$x, y = sm2$y)

      # 1st panel: Residuals vs Fitted
      p1 <- ggplot2::ggplot(df, ggplot2::aes(fitted, resid)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_line(data = sm1df, ggplot2::aes(x, y), inherit.aes = FALSE, colour = "red") +
        ggplot2::labs(
          title = "Residuals vs Fitted",
          subtitle = paste("Model:", deparse(.self$formula)),
          x = "Fitted values",
          y = "Residuals"
        ) +
        liu_theme()

      # 2nd panel: Scale Location
      p2 <- ggplot2::ggplot(df, ggplot2::aes(fitted, sl_y)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(data = sm2df, ggplot2::aes(x, y), inherit.aes = FALSE, colour = "red") +
        ggplot2::labs(
          title = "Scale\u2013Location",
          subtitle = paste("Model:", deparse(.self$formula)),
          x = "Fitted values",
          y = expression(sqrt("|Standardized residuals|"))
        ) +
        liu_theme()

      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required. Add it to Imports and install it.", call. = FALSE)
      }
      out <- p1 / p2

      base::print(out)
      invisible(out)
    },

    summary = function(...) {
      cat(sprintf(
        "linreg(formula = %s, data = %s)\n\n",
        deparse(.self$formula),
        .self$data_name
      ))
      cat("Coefficients:\n")
      printCoefmat(
        as.matrix(.self$summary_table),
        P.values = TRUE,
        has.Pvalue = TRUE
      )
      cat(
        sprintf(
          "Residual standard error: %.4f on %d degrees of freedom\n",
          .self$sigma,
          .self$df_residual
        )
      )
      invisible(.self)
    }
  )
)
