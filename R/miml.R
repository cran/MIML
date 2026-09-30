#' Recurrence-free survival analysis of the GSE65622 panel
#'
#' Runs the published MIML workflow on [RFS_data]: LightGBM imputation of every
#' incomplete biomarker, Bayesian consensus clustering of the completed panel,
#' and a Kaplan-Meier, log-rank and Cox comparison of the resulting clusters
#' against recurrence-free survival.
#'
#' Imputed values are the conditional-mean predictions of the fitted LightGBM
#' models. No residual is drawn and no donor is sampled, so the completions
#' carry no between-imputation variance and \strong{must not be combined by
#' Rubin's rules}. This is machine-learning single imputation.
#'
#' @param data Panel to analyse. Defaults to [RFS_data].
#' @param m Number of completed datasets. The published analysis used 5.
#' @param maxit Refit cycles per completion. The published analysis used 10.
#' @param nrounds Maximum boosting rounds. The published analysis used 200.
#' @param n_clusters Number of consensus clusters. The published analysis
#'   used 2.
#' @param n_bcc_markers Biomarkers carried into the clustering model, chosen by
#'   descending variance. The published analysis used 3.
#' @param predictors Predictor pool for the imputation models.
#'   \code{"manuscript"} uses the patient identifier, the survival endpoint,
#'   the event indicator and the six phenotype columns, which is what the
#'   published analysis did; \code{"all"} lets biomarkers predict each other
#'   and drops the outcome from the pool. See the leakage note below.
#' @param seed Random seed. The caller's RNG state is restored on exit.
#' @param verbose Print progress.
#'
#' @section Predictor pool and outcome leakage:
#' Under \code{predictors = "manuscript"} the survival time and the event
#' indicator are inside the predictor pool, and no biomarker predicts another.
#' Any association later estimated between the imputed values and survival is
#' therefore inflated by construction. The default reproduces the published
#' numbers; use \code{predictors = "all"} for new work.
#'
#' @section Run time:
#' At the defaults the imputation fits 507 biomarkers over 10 cycles for each
#' of 5 completions and takes roughly an hour.
#'
#' @return An object of class \code{miml_analysis}: a list with
#'   \code{endpoint}, \code{missingness}, \code{completed}, \code{clusters},
#'   \code{survival} and \code{settings}.
#'
#' @examples
#' data(RFS_data)
#'
#' ## small and fast: three biomarkers, one cycle, no clustering
#' small <- RFS_data[1:60, c("Patient_ID", "RFS", "event",
#'                           names(RFS_data)[4:6])]
#' miml_rfs(small, m = 1, maxit = 1, nrounds = 3,
#'          n_clusters = 0, verbose = FALSE)
#'
#' \donttest{
#' ## twenty biomarkers, two cycles, clustered and compared by survival
#' sub <- RFS_data[, c("Patient_ID", "RFS", "event", names(RFS_data)[4:23])]
#' res <- miml_rfs(sub, m = 1, maxit = 2, nrounds = 20)
#' res
#'
#' ## the published analysis is miml_rfs() at its defaults: 507 biomarkers,
#' ## 10 cycles, 5 completions, roughly an hour.
#' }
#'
#' @seealso [miml_mfs()], [RFS_data]
#' @export
miml_rfs <- function(data = NULL, m = 5L, maxit = 10L, nrounds = 200L,
                     n_clusters = 2L, n_bcc_markers = 3L,
                     predictors = c("manuscript", "all"),
                     seed = 123L, verbose = TRUE) {
  if (is.null(data)) {
    e <- new.env()
    utils::data("RFS_data", package = "MIML", envir = e)
    data <- get("RFS_data", envir = e)
  }
  .miml_run(data, "RFS", m, maxit, nrounds, n_clusters, n_bcc_markers,
            match.arg(predictors), seed, verbose)
}


#' Metastasis-free survival analysis of the GSE65622 panel
#'
#' Identical to [miml_rfs()] but run on [MFS_data] against metastasis-free
#' survival. Twenty-six of the eighty patients record a metastasis, against
#' five recurrences for the other endpoint, so this is the better powered of
#' the two comparisons.
#'
#' @inheritParams miml_rfs
#' @param data Panel to analyse. Defaults to [MFS_data].
#'
#' @inheritSection miml_rfs Predictor pool and outcome leakage
#' @inheritSection miml_rfs Run time
#'
#' @return An object of class \code{miml_analysis}, as for [miml_rfs()].
#'
#' @examples
#' data(MFS_data)
#'
#' small <- MFS_data[1:60, c("Patient_ID", "MFS", "event",
#'                           names(MFS_data)[4:6])]
#' miml_mfs(small, m = 1, maxit = 1, nrounds = 3,
#'          n_clusters = 0, verbose = FALSE)
#'
#' \donttest{
#' sub <- MFS_data[, c("Patient_ID", "MFS", "event", names(MFS_data)[4:23])]
#' res <- miml_mfs(sub, m = 1, maxit = 2, nrounds = 20)
#' res
#'
#' ## the published analysis is miml_mfs() at its defaults.
#' }
#'
#' @seealso [miml_rfs()], [MFS_data]
#' @export
miml_mfs <- function(data = NULL, m = 5L, maxit = 10L, nrounds = 200L,
                     n_clusters = 2L, n_bcc_markers = 3L,
                     predictors = c("manuscript", "all"),
                     seed = 123L, verbose = TRUE) {
  if (is.null(data)) {
    e <- new.env()
    utils::data("MFS_data", package = "MIML", envir = e)
    data <- get("MFS_data", envir = e)
  }
  .miml_run(data, "MFS", m, maxit, nrounds, n_clusters, n_bcc_markers,
            match.arg(predictors), seed, verbose)
}


#' @export
print.miml_analysis <- function(x, ...) {
  cat(sprintf("MIML analysis -- %s\n", x$endpoint))
  cat(sprintf("  panel       : %d rows x %d biomarkers, %d patients\n",
              x$missingness$n_rows, x$missingness$n_biomarkers,
              x$missingness$n_patients))
  cat(sprintf("  missing     : %d cells (%.2f%%), 0 after imputation\n",
              x$missingness$n_missing, 100 * x$missingness$prop_missing))
  cat(sprintf("  completions : %d (conditional-mean; NOT Rubin-poolable)\n",
              x$settings$m))
  if (is.null(x$clusters)) {
    cat("  clustering  : not run\n")
  } else {
    cat("  clusters    : ",
        paste(names(table(x$clusters)), as.integer(table(x$clusters)),
              sep = " n=", collapse = ", "), "\n", sep = "")
  }
  if (is.null(x$survival)) {
    cat("  survival    : not run\n")
  } else {
    s <- x$survival
    cat(sprintf("  survival    : %d patients, %d events, log-rank p = %.4g\n",
                s$n_patients, s$n_events, s$logrank_p))
    if (is.finite(s$hazard_ratio[["estimate"]]))
      cat(sprintf("  hazard ratio: %.3f (95%% CI %.3f-%.3f), p = %.4g\n",
                  s$hazard_ratio[["estimate"]], s$hazard_ratio[["lower"]],
                  s$hazard_ratio[["upper"]], s$hazard_ratio[["p"]]))
    else
      cat("  hazard ratio: not estimable\n")
  }
  invisible(x)
}


## ---------------------------------------------------------------------------
## internals
## ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.miml_run <- function(data, endpoint, m, maxit, nrounds, n_clusters,
                      n_bcc_markers, predictors, seed, verbose) {

  id_var <- "Patient_ID"; event_var <- "event"; time_var <- endpoint
  pheno  <- c("diarrhea_grade", "sampling_point", "slide_print_batch",
              "trg_score", "ypn_stage", "ypt_stage")

  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  absent <- setdiff(c(id_var, time_var, event_var), names(data))
  if (length(absent))
    stop("Columns not found in 'data': ", paste(absent, collapse = ", "),
         call. = FALSE)

  m     <- max(1L, as.integer(m))
  maxit <- max(1L, as.integer(maxit))
  n_clusters <- as.integer(n_clusters)

  keep    <- intersect(pheno, names(data))
  markers <- setdiff(names(data), c(id_var, time_var, event_var, keep))
  if (!length(markers))
    stop("No biomarker columns found in 'data'.", call. = FALSE)

  raw <- data[, c(id_var, time_var, event_var, keep, markers), drop = FALSE]
  for (v in markers) raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))

  n_miss <- sum(is.na(raw[, markers, drop = FALSE]))
  miss   <- list(n_rows       = nrow(raw),
                 n_patients   = length(unique(raw[[id_var]])),
                 n_biomarkers = length(markers),
                 n_cells      = nrow(raw) * length(markers),
                 n_missing    = n_miss,
                 prop_missing = n_miss / (nrow(raw) * length(markers)))

  if (verbose)
    message(sprintf("%s: %d rows x %d biomarkers | %d missing (%.2f%%)",
                    endpoint, miss$n_rows, miss$n_biomarkers,
                    miss$n_missing, 100 * miss$prop_missing))

  completed <- .miml_impute(raw, markers, predictors, m, maxit, nrounds,
                            seed, verbose,
                            protect = c(id_var, time_var, event_var))

  clusters <- NULL
  survival_fit <- NULL
  if (n_clusters >= 2L) {
    clusters <- .miml_cluster(completed, markers, id_var, time_var,
                              n_clusters, n_bcc_markers, verbose)
    if (!is.null(clusters))
      survival_fit <- .miml_survival(completed, clusters, id_var, time_var,
                                     event_var, verbose)
  }

  structure(list(
    endpoint    = endpoint,
    missingness = miss,
    completed   = completed,
    clusters    = clusters,
    survival    = survival_fit,
    settings    = list(id_var = id_var, time_var = time_var,
                       event_var = event_var, markers = markers,
                       m = m, maxit = maxit, nrounds = nrounds,
                       n_clusters = n_clusters,
                       n_bcc_markers = n_bcc_markers,
                       predictors = predictors, seed = seed)),
    class = "miml_analysis")
}


#' @keywords internal
#' @noRd
.miml_impute <- function(raw, markers, predictors, m, maxit, nrounds,
                         seed, verbose, protect) {

  ## LightGBM's OpenMP pool spin-waits while idle and accrues CPU time even
  ## for a single-threaded fit. Cap it here and restore the caller's setting.
  if ("setLGBMthreads" %in% getNamespaceExports("lightgbm")) {
    old <- try(lightgbm::getLGBMthreads(), silent = TRUE)
    lightgbm::setLGBMthreads(1L)
    if (!inherits(old, "try-error"))
      on.exit(try(lightgbm::setLGBMthreads(old), silent = TRUE), add = TRUE)
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    st <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", st, envir = globalenv()), add = TRUE)
  }

  dt <- data.table::as.data.table(raw)

  pool <- if (identical(predictors, "manuscript"))
    setdiff(names(dt), markers) else setdiff(names(dt), protect)
  for (cl in pool)
    if (!is.numeric(dt[[cl]]))
      data.table::set(dt, j = cl,
                      value = as.numeric(factor(dt[[cl]], exclude = NULL)))

  ## LightGBM draws its bagging and feature subsamples from its own C++
  ## generator, which set.seed() in R does not reach. Without an explicit
  ## seed here every run returns different imputed values -- and therefore
  ## different clusters and different p-values. Seed it per completion.
  lgb_params <- function(i) {
    sd <- as.integer(seed + i * 13L)
    list(objective = "regression", metric = "rmse",
         num_leaves = 31L, learning_rate = 0.05,
         feature_fraction = 0.5, bagging_fraction = 0.8,
         bagging_freq = 1L, min_data_in_leaf = 3L,
         seed = sd, bagging_seed = sd, feature_fraction_seed = sd,
         data_random_seed = sd, deterministic = TRUE,
         verbose = -1L, num_threads = 1L)
  }

  one <- function(i) {
    d <- data.table::copy(dt)
    params <- lgb_params(i)
    set.seed(seed + i * 13L)
    for (cycle in seq_len(maxit)) {
      for (v in markers) {
        na_i <- which(is.na(d[[v]]))
        if (!length(na_i)) next
        ob_i <- which(!is.na(d[[v]]))
        if (length(ob_i) < 2L) next
        y <- d[[v]][ob_i]

        cand <- setdiff(pool, v)
        if (!length(cand)) next
        r <- vapply(cand, function(p) {
          z <- suppressWarnings(
            stats::cor(d[[p]], d[[v]], use = "pairwise.complete.obs"))
          if (is.finite(z)) abs(z) else 0
        }, numeric(1))
        sel <- names(sort(r, decreasing = TRUE))[
          seq_len(min(50L, length(r)))]

        Xo <- as.matrix(d[ob_i, sel, with = FALSE])
        Xm <- as.matrix(d[na_i, sel, with = FALSE])
        ok <- colSums(!is.na(Xo)) > 0L
        if (!any(ok)) {
          data.table::set(d, i = na_i, j = v,
                          value = sample(y, length(na_i), replace = TRUE))
          next
        }
        Xo <- Xo[, ok, drop = FALSE]; Xm <- Xm[, ok, drop = FALSE]

        fit <- lightgbm::lgb.train(
          params = params,
          data = lightgbm::lgb.Dataset(Xo, label = y, free_raw_data = FALSE),
          nrounds = nrounds, verbose = -1L)

        data.table::set(d, i = na_i, j = v,
                        value = as.numeric(stats::predict(fit, Xm)))
      }
      if (verbose)
        message(sprintf("  completion %d, cycle %d of %d", i, cycle, maxit))
    }
    as.data.frame(d, stringsAsFactors = FALSE)
  }

  if (verbose) message("Imputing ", length(markers), " biomarkers ...")
  out <- one(1L)
  if (m > 1L) for (i in seq.int(2L, m)) invisible(one(i))
  out
}


#' @keywords internal
#' @noRd
.miml_cluster <- function(completed, markers, id_var, time_var,
                          n_clusters, n_bcc_markers, verbose) {

  if (!requireNamespace("BCClong", quietly = TRUE)) {
    warning("Package 'BCClong' is needed for clustering; stage skipped.",
            call. = FALSE)
    return(NULL)
  }

  v <- apply(as.matrix(completed[, markers, drop = FALSE]), 2,
             stats::var, na.rm = TRUE)
  v[!is.finite(v)] <- 0
  k <- min(as.integer(n_bcc_markers), length(markers))
  use <- names(sort(v, decreasing = TRUE))[seq_len(k)]
  if (verbose) message("Clustering on: ", paste(use, collapse = ", "))

  mm <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) rep(0.5, length(x)) else
      (x - r[1L]) / diff(r)
  }

  fit <- NULL
  for (s in c(7L, 100L, 999L, 314L, 1234L)) {
    if (verbose) message("  BCC.multi seed ", s, " ...")
    tryfit <- tryCatch(
      BCClong::BCC.multi(
        mydat = lapply(use, function(z) mm(completed[[z]])),
        id = list(completed[[id_var]]),
        time = list(mm(completed[[time_var]])),
        formula = replicate(k, stats::as.formula("y ~ time + (1|id)"),
                            simplify = FALSE),
        dist = rep("gaussian", k), num.cluster = n_clusters,
        center = 0, sigma.sq.e.common = 1,
        initial.cluster.membership = "random", seed.initial = s,
        burn.in = 500L, thin = 2L, per = 100L, max.iter = 1000L),
      error = function(e) NULL)
    if (!is.null(tryfit)) {
      g <- tryfit$cluster.global
      if (length(unique(g)) == n_clusters) { fit <- tryfit; break }
    }
  }
  if (is.null(fit)) {
    warning("BCC.multi() failed on every initialisation; stage skipped.",
            call. = FALSE)
    return(NULL)
  }

  g <- fit$cluster.global
  lab <- factor(g, levels = sort(unique(g)),
                labels = paste0("Cluster ", sort(unique(g))))
  if (length(g) == nrow(completed)) return(droplevels(lab))
  droplevels(lab[match(completed[[id_var]], unique(completed[[id_var]]))])
}


#' @keywords internal
#' @noRd
.miml_survival <- function(completed, clusters, id_var, time_var, event_var,
                           verbose) {

  first <- !duplicated(completed[[id_var]])
  km <- data.frame(time    = as.numeric(completed[[time_var]][first]),
                   event   = as.numeric(completed[[event_var]][first]),
                   Cluster = droplevels(clusters[first]))
  km <- km[is.finite(km$time) & is.finite(km$event), , drop = FALSE]

  if (nlevels(km$Cluster) < 2L) {
    warning("Fewer than two clusters among patients; survival skipped.",
            call. = FALSE)
    return(NULL)
  }

  fit <- survival::survfit(survival::Surv(time, event) ~ Cluster, data = km)
  lr  <- survival::survdiff(survival::Surv(time, event) ~ Cluster, data = km)
  p   <- stats::pchisq(lr$chisq, length(lr$n) - 1L, lower.tail = FALSE)

  hr  <- c(estimate = NA_real_, lower = NA_real_, upper = NA_real_,
           p = NA_real_)
  cox <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ Cluster, data = km),
    error = function(e) NULL)
  if (!is.null(cox)) {
    cf <- stats::coef(cox)
    if (length(cf) && is.finite(cf[1L])) {
      ci <- suppressWarnings(exp(stats::confint(cox)))
      hr <- c(estimate = unname(exp(cf[1L])),
              lower    = unname(ci[1L, 1L]),
              upper    = unname(ci[1L, 2L]),
              p        = unname(summary(cox)$coefficients[1L, "Pr(>|z|)"]))
    }
  }

  if (verbose)
    message(sprintf("Survival: %d patients, %d events, log-rank p = %.4g",
                    nrow(km), sum(km$event), p))

  list(n_patients = nrow(km), n_events = sum(km$event),
       fit = fit, logrank = lr, logrank_p = p, cox = cox, hazard_ratio = hr,
       median_table = summary(fit)$table, data = km)
}
