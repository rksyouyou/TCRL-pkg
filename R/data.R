#' A list containing TCR data, continuous reponse and covariate variables.
#'
#' A simulated dataset list containing the TCR information, covariate and response variables of 20 patients.
#'
#' @format A dataset containing the TCR information, covariate and response variables of 20 patients.
#' \describe{
#'   \item{TCRdat}{The TCR information with three columns, Subject.ID, AAseq, and Abundance, indicating the patient ID, amino acid sequence, and corresponding abundance.}
#'   \item{Y}{A vector of continous response.}
#'   \item{X}{A matrix of covariate variables.}
#'   \item{W}{A vector of Kyte and Doolittle hydrophobicity information.}
#' }
"example.data"
