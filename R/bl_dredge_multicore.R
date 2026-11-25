#' bl_dredge_multicore
#'
#' @description function built on MuMIn::dredge that makes maximum use of CPU power for faster dredging.
#'
#' @param global_model A global model (glm, glmm, lm, lmm, or any model type that MuMIn::dredge can handle)
#' @param packages packages that were used for fitting the global model. "lme4" is the standard here, but e.g. "lmerTest", "glmmTMB", "nlme", "MASS" (even combinations) are other candidates.
#' @param dredge_m.lim Vector of length == 2 specifying the min and max number of variables in each model. (see ?MuMIn::dredge --> m.lim). "Auto" sets to 1/10 of all cases of the global model.
#' @param dredge_subset equivalent to MuMIn::dredge --> subset
#' @param additional_dredge_arguments list of additional arguments accepted by MuMIn::dredge
#' @param free_cores number of cores that should not be used for parallel processing. Advice: leave at least one core out (for system performance). Modern ultra cores often have LP-E-cores (e.g.: ultra 5 235U has 2) that should be left out.

#' @export
#' @examples
#' # fit a global model, then do fast dredging and subsequent model averaging
#' the_independent_vars <- c("elevation", "summer_temperature", "winter_temperature", "summer_precipitation", "winter_precipitation", "snow_cover_duration")
#' the_global_model <- lme4::glmer(data = bl_ecology_data,
#'                                 formula = reformulate(c(the_independent_vars, "(1|sub_area)"), "grasshopper"),
#'                                 family = "poisson"
#' )
#' summary(the_global_model)
#' grasshopper_dredge <- bl_dredge_multicore(the_global_model)
#' grasshopper_dredge
#' MuMIn::model.avg(grasshopper_dredge,  subset = delta < 1)



bl_dredge_multicore <- function(global_model,
                                packages = c("lme4"),
                                dredge_m.lim = "auto",
                                dredge_subset = NULL,
                                additional_dredge_arguments = list(),
                                free_cores = 2) {

    envir <- parent.frame()

    # Subset-Ausdruck unausgewertet erfassen (VOR jedem Zugriff auf dredge_subset!)
    subset_expr <- substitute(dredge_subset)
    has_subset  <- !identical(subset_expr, quote(NULL))

    old_opt <- options(na.action = "na.fail")
    on.exit(options(old_opt), add = TRUE)

    if (identical(dredge_m.lim, "auto")) {
        dredge_m.lim <- c(0, floor(nobs(global_model) / 10))
        message("dredge_m.lim automatisch gesetzt: c(", dredge_m.lim[1], ", ", dredge_m.lim[2], ")")
    }

    dredge_args <- c(
        list(global_model, evaluate = FALSE),
        if (!is.null(dredge_m.lim)) list(m.lim = dredge_m.lim),
        if (has_subset)             list(subset = subset_expr),
        additional_dredge_arguments
    )
    calls <- do.call(MuMIn::dredge, dredge_args)
    # envir <- parent.frame()
    #
    # old_opt <- options(na.action = "na.fail")
    # on.exit(options(old_opt), add = TRUE)
    #
    # # --- 1. Modellaufrufe erzeugen (fittet nichts) ---------------------
    # # "auto": maximal 1/10 der Beobachtungen, NULL: kein Limit, sonst Vektor c(min, max)
    # if (identical(dredge_m.lim, "auto")) {
    #     dredge_m.lim <- c(0, floor(nobs(global_model) / 10))
    #     message("dredge_m.lim automatisch gesetzt: c(", dredge_m.lim[1], ", ", dredge_m.lim[2], ")")
    # }
    #
    # dredge_args <- c(
    #     list(global_model, evaluate = FALSE),
    #     if (!is.null(dredge_m.lim))  list(m.lim  = dredge_m.lim),
    #     if (!is.null(dredge_subset)) list(subset = dredge_subset),
    #     additional_dredge_arguments
    # )
    # calls <- do.call(MuMIn::dredge, dredge_args)
    message(length(calls), " Modelle werden gefittet")

    # --- 2. Cluster ------------------------------------------------------
    n_cores <- max(1, parallel::detectCores() - free_cores)
    cl <- parallel::makeCluster(n_cores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    packages <- c("MuMIn", packages)
    parallel::clusterExport(cl, "packages", envir = environment())
    parallel::clusterEvalQ(cl, lapply(packages, library, character.only = TRUE))

    # Alle Objekte, die im Call des globalen Modells vorkommen, exportieren
    needed <- all.vars(getCall(global_model))
    needed <- needed[vapply(needed, exists, logical(1), envir = envir)]
    parallel::clusterExport(cl, needed, envir = envir)

    fit_fun <- function(x) {
        tryCatch(eval(x, envir = globalenv()),
                 error = function(e) conditionMessage(e))
    }
    environment(fit_fun) <- globalenv()

    # --- 3. Fitten --------------------------------------------------------
    fits <- parallel::parLapplyLB(cl, calls, fit_fun#, chunk.size = max(1, ceiling(length(calls) / (n_cores * 20)))
                                  )
    parallel::stopCluster(cl)

    # --- 4. Fehlgeschlagene Fits -----------------------------------------
    ok <- !vapply(fits, is.character, logical(1))
    if (any(!ok)) {
        warning(sum(!ok), " Fits fehlgeschlagen: ",
                paste(unique(unlist(fits[!ok])), collapse = " | "))
    }
    failed_calls <- calls[!ok]
    fits <- fits[ok]

    dredge_object <- MuMIn::model.sel(fits)

    # --- 5. "Normales" dredge-Objekt nachbauen ----------------------------
    attr(dredge_object, "modelList")    <- NULL
    attr(dredge_object, "model.family") <- NULL
    attr(dredge_object, "global")       <- global_model
    attr(dredge_object, "global.call")  <- getCall(global_model)
    attr(dredge_object, "rank")         <- MuMIn::AICc
    attr(dredge_object, "failed_calls") <- failed_calls
    dredge_object
}

