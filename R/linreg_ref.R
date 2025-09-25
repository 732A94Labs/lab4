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
      # The test "print() method works" strips context, so match.call() ("data_label <- deparse(call[[3]]) ") fails
      # to find the original variable name (for the test: "iris") of the dataframe. This is a kind of sloppy workaround that performs
      # a "reverse lookup" to manually search the parent environments to find the name
      # of the variable that is identical to the provided data object.
      data_label <- ""
      for (i in 1:sys.nframe()) {
        env <- parent.frame(n = i)
        vars <- ls(env)
        for (var in vars) {
          if (identical(get(var, env), data)) {
            data_label <- var
            break
          }
        }
        if (data_label != "") break
      }

      mf <- model.frame(formula, data)
      y <- model.response(mf)
      X <- model.matrix(attr(mf, "terms"), mf)

      # beta_hat = (X'X)^(-1) X'y (regression coefficients)
      coefficients <- as.numeric(solve(crossprod(X)) %*% crossprod(X, y))
      names(coefficients) <- colnames(X)

      # y_hat = X * beta_hat (fitted values)
      fitted <- as.numeric(X %*% coefficients)

      # residuals = y - y_hat
      residuals <- as.numeric(y - fitted)

      # df = n - p (degrees of freedom)
      df <- length(y) - ncol(X)
      if (df <= 0) {
        stop("Residual degrees of freedom must be positive.")
      }

      # sigma_hat^2 = RSS / df (residual standard error)
      rss <- crossprod(residuals)
      sigma_hat_squared <- as.numeric(rss / df)

      # var of resids = sigma_hat^2 * (X'X)^(-1) (covariance matrix of coefficients)'
      var_hat_beta <- sigma_hat_squared * solve(crossprod(X))

      # t_beta = beta_hat / se(beta_hat) (t-values and p-values)
      t_values <- coefficients / sqrt(diag(var_hat_beta))
      p_values <- 2 * pt(-abs(t_values), df)

      coef_table <- data.frame(
        Estimate = coefficients,
        `Std. Error` = sqrt(diag(var_hat_beta)),
        `t value` = t_values,
        `Pr(>|t|)` = p_values,
        check.names = FALSE
      )

      .self$formula <- formula
      .self$data <- data
      .self$data_name <- data_label
      .self$coefficients <- coefficients
      .self$fitted_values <- fitted
      .self$residuals <- residuals
      .self$df_residual <- df
      .self$sigma <- sqrt(sigma_hat_squared)
      .self$var_hat_beta <- var_hat_beta
      .self$summary_table <- coef_table

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
