#' Recurrence-free survival panel (GSE65622)
#'
#' Longitudinal serum proteomic measurements from patients with locally
#' advanced rectal cancer undergoing multimodal treatment, paired with
#' recurrence-free survival. The default input to [miml_rfs()].
#'
#' @format A data frame of 288 rows and 516 columns: 80 patients observed at
#'   three or four sampling points each. \code{Patient_ID} is an integer
#'   identifier repeated across a patient's visits; \code{RFS} is
#'   recurrence-free survival in days measured \emph{from that visit}, so it
#'   decreases across a patient's rows; \code{event} is the recurrence
#'   indicator, constant within a patient. Then 507 protein biomarkers named
#'   as in the original submission, then six phenotype columns:
#'   \code{diarrhea_grade}, \code{sampling_point}, \code{slide_print_batch},
#'   \code{trg_score}, \code{ypn_stage}, \code{ypt_stage}.
#'
#'   Biomarker names are not syntactic R names, so quote them or use
#'   \code{[[ ]]}. About 10.3 percent of the biomarker measurements are
#'   missing. Five of the 80 patients record a recurrence.
#'
#' @source Gene Expression Omnibus accession
#'   \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65622}{GSE65622},
#'   released April 2016.
#'
#' @seealso [MFS_data], [miml_rfs()]
#'
#' @examples
#' data(RFS_data)
#' dim(RFS_data)
#' RFS_data[1:5, 1:6]
"RFS_data"


#' Metastasis-free survival panel (GSE65622)
#'
#' The same patients, visits and 507 biomarkers as [RFS_data], paired instead
#' with metastasis-free survival. The default input to [miml_mfs()].
#'
#' @format A data frame of 288 rows and 516 columns, identical in structure to
#'   [RFS_data] except that the survival column is \code{MFS} and \code{event}
#'   is the metastasis indicator. Twenty-six of the 80 patients record a
#'   metastasis.
#'
#' @source Gene Expression Omnibus accession
#'   \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65622}{GSE65622},
#'   released April 2016.
#'
#' @seealso [RFS_data], [miml_mfs()]
#'
#' @examples
#' data(MFS_data)
#' dim(MFS_data)
#' MFS_data[1:5, 1:6]
"MFS_data"
