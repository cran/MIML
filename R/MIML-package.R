#' MIML: machine-learning imputation, clustering and survival for GSE65622
#'
#' Two functions run the published MIML workflow end to end: [miml_rfs()] on
#' [RFS_data] against recurrence-free survival, and [miml_mfs()] on [MFS_data]
#' against metastasis-free survival. Each performs LightGBM imputation of the
#' incomplete biomarker panel, Bayesian consensus clustering of the completed
#' panel, and a survival comparison of the resulting clusters.
#'
#' Imputed values are conditional-mean predictions, so the completions carry
#' no between-imputation variance and \strong{must not be combined by Rubin's
#' rules}. This is machine-learning single imputation.
#'
#' @references
#' Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., Ye, Q. and
#' Liu, T.-Y. (2017). LightGBM: a highly efficient gradient boosting decision
#' tree. \emph{Advances in Neural Information Processing Systems} \strong{30},
#' 3146-3154.
#' \url{https://papers.nips.cc/paper/6907-lightgbm-a-highly-efficient-gradient-boosting-decision-tree}
#'
#' Lock, E. F. and Dunson, D. B. (2013). Bayesian consensus clustering.
#' \emph{Bioinformatics} \strong{29}(20), 2610-2616.
#' \doi{10.1093/bioinformatics/btt425}
#'
#' @keywords internal
"_PACKAGE"
