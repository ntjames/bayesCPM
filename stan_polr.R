function (formula, data, weights, ..., subset, na.action = getOption("na.action",
    "na.omit"), contrasts = NULL, model = TRUE, method = c("logistic",
    "probit", "loglog", "cloglog", "cauchit"), prior = R2(stop("'location' must be specified")),
    prior_counts = dirichlet(1), shape = NULL, rate = NULL, prior_PD = FALSE,
    algorithm = c("sampling", "meanfield", "fullrank"), adapt_delta = NULL,
    do_residuals = NULL)
{
    data <- validate_data(data, if_missing = environment(formula))
    algorithm <- match.arg(algorithm)
    if (is.null(do_residuals))
        do_residuals <- algorithm == "sampling"
    call <- match.call(expand.dots = TRUE)
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    if (is.matrix(eval.parent(m$data))) {
        m$data <- as.data.frame(data)
    }
    else {
        m$data <- data
    }
    m$method <- m$model <- m$... <- m$prior <- m$prior_counts <- m$prior_PD <- m$algorithm <- m$adapt_delta <- m$shape <- m$rate <- m$do_residuals <- NULL
    m[[1L]] <- quote(stats::model.frame)
    m$drop.unused.levels <- FALSE
    m <- eval.parent(m)
    m <- check_constant_vars(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0L) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1L
    }
    else stop("an intercept is needed and assumed")
    K <- ncol(x)
    wt <- model.weights(m)
    if (!length(wt))
        wt <- rep(1, n)
    offset <- model.offset(m)
    if (length(offset) <= 1L)
        offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y))
        stop("Response variable must be a factor.", call. = FALSE)
    lev <- levels(y)
    llev <- length(lev)
    if (llev < 2L)
        stop("Response variable must have 2 or more levels.",
            call. = FALSE)
    q <- llev - 1L
    stanfit <- stan_polr.fit(x = x, y = y, wt = wt, offset = offset,
        method = method, prior = prior, prior_counts = prior_counts,
        shape = shape, rate = rate, prior_PD = prior_PD, algorithm = algorithm,
        adapt_delta = adapt_delta, do_residuals = do_residuals,
        ...)
    if (algorithm != "optimizing" && !is(stanfit, "stanfit"))
        return(stanfit)
    inverse_link <- linkinv(method)
    if (llev == 2L) {
        family <- switch(method, logistic = binomial(link = "logit"),
            loglog = binomial(loglog), binomial(link = method))
        fit <- nlist(stanfit, family, formula, offset, weights = wt,
            x = cbind(`(Intercept)` = 1, x), y = as.integer(y ==
                lev[2]), data, call, terms = Terms, model = m,
            algorithm, na.action = attr(m, "na.action"), contrasts = attr(x,
                "contrasts"), stan_function = "stan_polr")
        out <- stanreg(fit)
        if (!model)
            out$model <- NULL
        if (algorithm == "sampling")
            check_rhats(out$stan_summary[, "Rhat"])
        if (is.null(shape) && is.null(rate))
            return(out)
        out$method <- method
        return(structure(out, class = c("stanreg", "polr")))
    }
    K2 <- K + llev - 1
    stanmat <- as.matrix(stanfit)[, 1:K2, drop = FALSE]
    covmat <- cov(stanmat)
    coefs <- apply(stanmat[, 1:K, drop = FALSE], 2L, median)
    ses <- apply(stanmat[, 1:K, drop = FALSE], 2L, mad)
    zeta <- apply(stanmat[, (K + 1):K2, drop = FALSE], 2L, median)
    eta <- linear_predictor(coefs, x, offset)
    mu <- inverse_link(eta)
    means <- rstan::get_posterior_mean(stanfit)
    residuals <- means[grep("^residuals", rownames(means)), ncol(means)]
    names(eta) <- names(mu) <- rownames(x)
    if (!prior_PD) {
        if (!do_residuals) {
            residuals <- rep(NA, times = n)
        }
        names(residuals) <- rownames(x)
    }
    stan_summary <- make_stan_summary(stanfit)
    if (algorithm == "sampling")
        check_rhats(stan_summary[, "Rhat"])
    out <- nlist(coefficients = coefs, ses, zeta, residuals,
        fitted.values = mu, linear.predictors = eta, covmat,
        y, x, model = if (model)
            m, data, offset, weights = wt, prior.weights = wt,
        family = method, method, contrasts, na.action, call,
        formula, terms = Terms, prior.info = attr(stanfit, "prior.info"),
        algorithm, stan_summary, stanfit, rstan_version = utils::packageVersion("rstan"),
        stan_function = "stan_polr")
    structure(out, class = c("stanreg", "polr"))
}
<bytecode: 0x5574a8fcdf98>
<environment: namespace:rstanarm>
