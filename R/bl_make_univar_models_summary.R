

#' bl_make_univar_models_summary
#'
#' @description bl function to make quick overview tables of a set of univariable
#' models (GLM or GLMM) using the same dependent variable and different independent
#' variables. All models use the same family. Models can be fitted with base `glm`,
#' `lme4` or `glmmTMB`. For every model, DHARMa residual diagnostics (dispersion,
#' zero-inflation, outliers) can be run; problems are reported in a column and as a warning.
#'
#' @param data_table a dataframe with a dependent and at least one independent variable
#' @param dependent_var a character string of the dependent variable name
#' @param independent_vars a vector of character strings of the independent variable names
#' @param random_formula optional: if a GLMM should be calculated, a character string
#'   in form of a random formula, e.g. "(1|site)"
#' @param offset_var optional: a character string of the offset variable name (will be log-transformed)
#' @param model_family a family. Can be a character string ("poisson", "nbinom2", "beta"),
#'   a family function (poisson, glmmTMB::nbinom2) or a family object (Gamma(link = "log"),
#'   glmmTMB::ordbeta()). Use an object if you need a non-default link.
#' @param family deprecated, use `model_family`
#' @param engine "auto" (default), "lme4" or "glmmTMB". "auto" uses glmmTMB for families only available there
#'   (nbinom1/2, beta, ordbeta, tweedie, compois, ...) or if `zi_formula`/`disp_formula` is given,
#'   otherwise lme4 (GLMM) or glm (no random formula). "glmmTMB" can also be used without random effects.
#' @param zi_formula optional, glmmTMB only: zero-inflation formula, e.g. `~1`
#' @param disp_formula optional, glmmTMB only: dispersion formula, e.g. `~1`
#' @param null_model for R² calculation, a null model can be specified. Especially useful
#'   for mixed models, where conditional R² computation may require a null model.
#' @param tc_off logical. if the family is binomial with 0/1 response, tc_off turns off
#'   separate columns for treatment and control cases
#' @param do_scale logical. should numerical independent variables be scaled for Estimate calculation?
#' @param P_star logical. should P values be translated into "**" stars?
#' @param mean logical. If true, variable means ± standard errors are added to the output table
#' @param error_type either "SE" or "SD"
#' @param R_2 logical. should r-squared values be added to the output table?
#' @param AUC logical. should AUC values be added to the output table (binomial only)?
#' @param diagnostics logical. run DHARMa residual diagnostics for every model?
#' @param diag_tests which tests to run: any of "dispersion", "zeroinflation", "outliers", "uniformity"
#' @param diag_alpha significance level for the diagnostic tests
#' @param diag_n_sim number of simulations for DHARMa::simulateResiduals
#' @param diag_thresholds list defining when a significant test counts as *severe*:
#'   `dispersion` = c(lower, upper) bounds of the dispersion ratio,
#'   `zeroinflation` = upper bound of the observed/expected zero ratio.
#'   Significant outlier excess is always treated as severe.
#' @param diag_column logical. add a "Diagnostics" column to the result table?
#'   Severe problems are marked with "[!]". The complete test statistics are always
#'   stored in `attr(result, "diagnostics")`.
#' @param diag_warn logical. emit one warning listing all models with severe problems?
#' @param debugging logical. print intermediate objects
#' @param max_digits integer. indicates the maximum of digits (up to 4).
#'
#' @examples
#' # lme4 (as before)
#' bl_make_univar_models_summary(bl_ecology_data, "grasshopper",
#'                               names(bl_ecology_data)[4:20],
#'                               "(1|sub_area)", model_family = "poisson")
#'
#' # glmmTMB with negative binomial
#' bl_make_univar_models_summary(bl_ecology_data, "grasshopper",
#'                               names(bl_ecology_data)[4:20],
#'                               "(1|sub_area)", model_family = "nbinom2",
#'                               engine = "glmmTMB")
#'
#' # glmmTMB, Gamma with log link, no random effect
#' bl_make_univar_models_summary(dat, "shannon", c("elevation", "area"),
#'                               model_family = Gamma(link = "log"),
#'                               engine = "glmmTMB")
#'
#' # zero-inflation model
#' bl_make_univar_models_summary(dat, "richness", c("elevation", "area"),
#'                               "(1|block)", model_family = "nbinom2",
#'                               zi_formula = ~1)
#'
#' @export
bl_make_univar_models_summary <-
    function(data_table,
             dependent_var,
             independent_vars,
             random_formula = NA,
             offset_var = NA,
             family = NULL,
             model_family = NULL,
             engine = c("auto", "lme4", "glmmTMB"),
             zi_formula = NULL,
             disp_formula = NULL,
             null_model = NULL,
             tc_off = FALSE,
             do_scale = TRUE,
             P_star = TRUE,
             mean = TRUE,
             error_type = c("SE", "SD"),
             R_2 = TRUE,
             AUC = FALSE,
             diagnostics = TRUE,
             diag_tests = c("dispersion", "zeroinflation", "outliers"),
             diag_alpha = 0.05,
             diag_n_sim = 250,
             diag_thresholds = list(dispersion = c(0.5, 1.5), zeroinflation = 1.5),
             diag_column = TRUE,
             diag_warn = TRUE,
             debugging = FALSE,
             max_digits = 3
    ) {

        error_type <- match.arg(error_type)
        engine <- match.arg(engine)
        diag_thresholds <- utils::modifyList(
            list(dispersion = c(0.5, 1.5), zeroinflation = 1.5), diag_thresholds)

        dbg <- function(label, obj = NULL) {
            if (isTRUE(debugging)) {
                cat("\n\033[31m", label, "\033[0m\n", sep = "")
                if (!is.null(obj)) print(obj)
            }
        }

        # ---- family ---------------------------------------------------------
        if (is.null(model_family) && !is.null(family)) {
            message("'family' is deprecated. Please use 'model_family' instead.")
            model_family <- family
        }
        if (is.null(model_family)) stop("model_family must be specified")
        fam <- .bl_resolve_family(model_family)
        dbg("family:", fam)

        # ---- random part / engine ---------------------------------------------
        has_random <- !is.null(random_formula) && !is.na(random_formula)
        if (has_random && !grepl("^\\(.*\\)$", random_formula)) {
            random_formula <- paste0("(", random_formula, ")")
        }

        uses_tmb_features <- !is.null(zi_formula) || !is.null(disp_formula)
        if (engine == "auto") {
            engine <- if (fam$family %in% .bl_tmb_only_families || uses_tmb_features) "glmmTMB" else "lme4"
        }
        if (engine == "lme4" && (fam$family %in% .bl_tmb_only_families || uses_tmb_features)) {
            stop("Family '", fam$family, "' / zi_formula / disp_formula require engine = 'glmmTMB'.")
        }
        if (engine == "glmmTMB" && !requireNamespace("glmmTMB", quietly = TRUE)) {
            stop("Package 'glmmTMB' is required for engine = 'glmmTMB'.")
        }
        if (engine == "lme4" && has_random && !requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for random effects with engine = 'lme4'.")
        }
        if (diagnostics && !requireNamespace("DHARMa", quietly = TRUE)) {
            warning("Package 'DHARMa' not installed: diagnostics are skipped.")
            diagnostics <- FALSE
        }
        message("Engine: ", if (engine == "lme4" && !has_random) "glm" else engine,
                " | family: ", fam$family, " (", fam$link, ")")

        # ---- data / offset ----------------------------------------------------
        dat <- as.data.frame(data_table)
        use_offset <- !is.null(offset_var) && !is.na(offset_var)
        if (use_offset) dat$.bl_offset <- log(dat[[offset_var]])

        make_formula <- function(fixed_part) {
            rhs <- c(fixed_part, if (has_random) random_formula, if (use_offset) "offset(.bl_offset)")
            stats::as.formula(paste(dependent_var, "~", paste(rhs, collapse = " + ")))
        }

        fit_model <- function(formula) {
            if (engine == "glmmTMB") {
                zi <- if (is.null(zi_formula)) ~0 else zi_formula
                dp <- if (is.null(disp_formula)) ~1 else disp_formula
                glmmTMB::glmmTMB(formula, data = dat, family = fam,
                                 ziformula = zi, dispformula = dp)
            } else if (has_random) {
                if (fam$family == "gaussian" && fam$link == "identity") {
                    lme4::lmer(formula, data = dat, REML = FALSE)
                } else {
                    lme4::glmer(formula, data = dat, family = fam)
                }
            } else {
                stats::glm(formula, data = dat, family = fam)
            }
        }

        # ---- reference fit ("0-model") -----------------------------------------
        ref_fit <- fit_model(make_formula("1"))
        dbg("reference fit (or '0-model'):", ref_fit)

        # ---- binary response? (treatment / control means) ----------------------
        is_bin <- fam$family == "binomial" &&
            dependent_var %in% names(dat) &&
            all(dat[[dependent_var]] %in% c(0, 1, NA)) &&
            !tc_off

        fmt_mean_err <- function(x) {
            x <- x[!is.na(x)]
            err <- if (error_type == "SE") stats::sd(x) / sqrt(length(x)) else stats::sd(x)
            paste0(bl_round_beautifully(base::mean(x), last = FALSE, max_dgts = max_digits),
                   " ± ",
                   bl_round_beautifully(err, last = TRUE, max_dgts = max_digits))
        }

        results <- data.frame("Variable" = independent_vars)
        diag_list <- list()
        severe_msgs <- character()

        # ---- loop over independent variables -----------------------------------
        for (my_var in independent_vars) {

            dbg("Variable in process:", my_var)
            row_i <- results$Variable == my_var
            is_num <- is.numeric(dat[[my_var]])

            ## i) means -----------------------------------------------------------
            if (isTRUE(mean)) {
                mean_col <- paste0("Mean ± ", error_type)
                if (is_bin) {
                    y <- dat[[dependent_var]]
                    tcol <- paste0("Treatment ", mean_col)
                    ccol <- paste0("Control ", mean_col)
                    if (is_num) {
                        results[row_i, tcol] <- fmt_mean_err(dat[[my_var]][which(y == 1)])
                        results[row_i, ccol] <- fmt_mean_err(dat[[my_var]][which(y == 0)])
                    } else {
                        results[row_i, tcol] <- "."
                        results[row_i, ccol] <- "."
                    }
                } else {
                    results[row_i, mean_col] <- if (is_num) fmt_mean_err(dat[[my_var]]) else "."
                }
            }

            ## ii) univariable model ------------------------------------------------
            fixed_term <- if (is_num && isTRUE(do_scale)) paste0("scale(", my_var, ")") else my_var
            my_formula <- make_formula(fixed_term)
            dbg("Formula:", my_formula)

            uv_fit <- fit_model(my_formula)
            dbg("Model summary:", summary(uv_fit))

            ## p value (likelihood ratio test against the null model)
            if (stats::nobs(ref_fit) != stats::nobs(uv_fit)) {
                warning("'", my_var, "': number of observations differs from the null model ",
                        "(missing values?). The LRT p value is not valid.", call. = FALSE)
            }
            an <- if (engine == "lme4" && !has_random) {
                stats::anova(ref_fit, uv_fit, test = "Chisq")
            } else {
                suppressMessages(stats::anova(ref_fit, uv_fit))
            }
            p_col <- grep("^Pr\\(", names(an), value = TRUE)
            p_value <- round(an[[p_col[1]]][2], 3)
            dbg("P value:", p_value)

            if (isTRUE(P_star)) {
                p_value <- if (is.na(p_value)) "n.a."
                else if (p_value < 0.001) "***"
                else if (p_value < 0.01) "**"
                else if (p_value < 0.05) "*"
                else "n.s."
            }

            ## estimate and SE (first non-intercept coefficient of the conditional model)
            cf <- .bl_coef_table(uv_fit)
            est <- bl_round_beautifully(cf[2, 1], last = FALSE)
            std_err <- bl_round_beautifully(cf[2, 2], last = TRUE)
            results[row_i, "Estimate ± SE"] <- paste0(est, " ± ", std_err)
            results[row_i, "P"] <- p_value

            ## R² -----------------------------------------------------------------
            if (isTRUE(R_2)) {
                my_R2 <- tryCatch(
                    if (!is.null(null_model)) performance::r2(uv_fit, null_model = null_model)
                    else performance::r2(uv_fit),
                    error = function(e) {
                        warning("R² could not be computed for '", my_var, "': ",
                                conditionMessage(e), call. = FALSE)
                        NULL
                    })
                dbg("R²:", my_R2)
                if (is.null(my_R2)) {
                    results[row_i, "R²"] <- NA
                } else if (length(my_R2) == 1) {
                    results[row_i, "R²"] <- round(my_R2[[1]], 2)
                } else if (length(my_R2) == 2) {
                    rc <- if (is.na(my_R2[[1]])) my_R2[[2]] else my_R2[[1]]
                    results[row_i, "Rc/Rm"] <- paste0(round(rc, 2), "/", round(my_R2[[2]], 2))
                }
            }

            ## AUC ----------------------------------------------------------------
            if (isTRUE(AUC)) {
                if (fam$family == "binomial") {
                    prob <- as.numeric(stats::predict(uv_fit, type = "response"))
                    obs <- insight::get_response(uv_fit)
                    roc_obj <- pROC::roc(response = obs, predictor = prob, quiet = TRUE)
                    results[row_i, "AUC"] <- as.numeric(pROC::auc(roc_obj))
                } else {
                    warning("AUC only implemented for binomial models. No AUC calculated")
                }
            }

            ## iii) residual diagnostics --------------------------------------------
            if (isTRUE(diagnostics)) {
                y_obs <- insight::get_response(uv_fit)
                dg <- .bl_diagnose(uv_fit,
                                   y_has_zero = any(y_obs == 0, na.rm = TRUE),
                                   tests = diag_tests, alpha = diag_alpha,
                                   n_sim = diag_n_sim, thr = diag_thresholds)
                diag_list[[my_var]] <- data.frame(Variable = my_var, dg$stats, flags = dg$text,
                                                  severe = length(dg$severe) > 0,
                                                  stringsAsFactors = FALSE)
                if (isTRUE(diag_column)) results[row_i, "Diagnostics"] <- dg$text
                if (length(dg$severe) > 0) {
                    severe_msgs <- c(severe_msgs,
                                     paste0("  - ", my_var, ": ", paste(dg$severe, collapse = ", ")))
                }
            }
        }

        if (isTRUE(diagnostics)) {
            attr(results, "diagnostics") <- do.call(rbind, diag_list)
            if (isTRUE(diag_warn) && length(severe_msgs) > 0) {
                warning("Severe model diagnostic problems (DHARMa, alpha = ", diag_alpha, "):\n",
                        paste(severe_msgs, collapse = "\n"),
                        "\nSee column 'Diagnostics' / attr(result, 'diagnostics'). ",
                        "Consider another family, a dispersion/zero-inflation formula, or checking outliers.",
                        call. = FALSE)
            }
        }

        dbg("Result table:")
        return(results)
    }


# ---- internal helpers ---------------------------------------------------------

# families that only exist in glmmTMB (value of fam$family)
.bl_tmb_only_families <- c("nbinom1", "nbinom2", "beta", "ordbeta", "tweedie", "compois",
                           "genpois", "betabinomial", "lognormal", "t", "skewnormal",
                           "truncated_poisson", "truncated_nbinom1", "truncated_nbinom2")

# exported glmmTMB family function names (for character input)
.bl_tmb_family_fns <- c("nbinom1", "nbinom2", "beta_family", "ordbeta", "tweedie", "compois",
                        "genpois", "betabinomial", "lognormal", "t_family", "skewnormal",
                        "truncated_poisson", "truncated_nbinom1", "truncated_nbinom2")

#' @noRd
.bl_resolve_family <- function(x) {
    if (inherits(x, "family")) return(x)
    if (is.function(x)) return(x())
    if (is.character(x) && length(x) == 1) {
        nm <- if (x == "beta") "beta_family" else x
        if (nm %in% .bl_tmb_family_fns) {
            if (!requireNamespace("glmmTMB", quietly = TRUE)) {
                stop("Family '", x, "' requires package 'glmmTMB'.")
            }
            return(getExportedValue("glmmTMB", nm)())
        }
        return(get(nm, mode = "function", envir = asNamespace("stats"))())
    }
    stop("model_family must be a character string, a family function or a family object.")
}

# coefficient matrix of the conditional model (glmmTMB returns a list)
#' @noRd
.bl_coef_table <- function(fit) {
    cf <- summary(fit)$coefficients
    if (is.list(cf) && !is.matrix(cf)) cf <- cf$cond
    cf
}

# DHARMa diagnostics for one model; returns test statistics, a flag text and severe problems
#' @noRd
.bl_diagnose <- function(fit, y_has_zero, tests, alpha, n_sim, thr) {

    stats_out <- data.frame(disp_ratio = NA_real_, disp_p = NA_real_,
                            zi_ratio = NA_real_, zi_p = NA_real_,
                            n_outliers = NA_real_, outlier_p = NA_real_,
                            uniformity_p = NA_real_)
    flags <- character()
    severe <- character()
    add_flag <- function(txt, is_severe) {
        flags <<- c(flags, paste0(txt, if (is_severe) " [!]" else ""))
        if (is_severe) severe <<- c(severe, txt)
    }

    res <- tryCatch(DHARMa::simulateResiduals(fit, n = n_sim, plot = FALSE),
                    error = function(e) NULL)
    if (is.null(res)) {
        return(list(stats = stats_out, text = "n/a (DHARMa failed)", severe = character()))
    }

    if ("dispersion" %in% tests) {
        t <- tryCatch(DHARMa::testDispersion(res, plot = FALSE), error = function(e) NULL)
        if (!is.null(t)) {
            r <- unname(t$statistic)
            stats_out$disp_ratio <- r
            stats_out$disp_p <- t$p.value
            if (!is.na(t$p.value) && t$p.value < alpha) {
                over <- r > 1
                add_flag(sprintf("%s (ratio %.2f)", if (over) "overdispersion" else "underdispersion", r),
                         if (over) r > thr$dispersion[2] else r < thr$dispersion[1])
            }
        }
    }

    if ("zeroinflation" %in% tests && y_has_zero) {
        t <- tryCatch(DHARMa::testZeroInflation(res, plot = FALSE), error = function(e) NULL)
        if (!is.null(t)) {
            r <- unname(t$statistic)
            stats_out$zi_ratio <- r
            stats_out$zi_p <- t$p.value
            if (!is.na(t$p.value) && t$p.value < alpha && r > 1) {
                add_flag(sprintf("zero-inflation (obs/exp %.2f)", r), r > thr$zeroinflation)
            }
        }
    }

    if ("outliers" %in% tests) {
        t <- tryCatch(DHARMa::testOutliers(res, alternative = "greater", plot = FALSE),
                      error = function(e) NULL)
        if (!is.null(t)) {
            stats_out$n_outliers <- unname(t$statistic)
            stats_out$outlier_p <- t$p.value
            if (!is.na(t$p.value) && t$p.value < alpha) {
                add_flag(sprintf("outliers (n = %d)", as.integer(t$statistic)), TRUE)
            }
        }
    }

    if ("uniformity" %in% tests) {
        t <- tryCatch(DHARMa::testUniformity(res, plot = FALSE), error = function(e) NULL)
        if (!is.null(t)) {
            stats_out$uniformity_p <- t$p.value
            if (!is.na(t$p.value) && t$p.value < alpha) {
                add_flag(sprintf("residual distribution (KS p = %.3f)", t$p.value), FALSE)
            }
        }
    }

    list(stats = stats_out,
         text = if (length(flags) == 0) "OK" else paste(flags, collapse = "; "),
         severe = severe)
}
