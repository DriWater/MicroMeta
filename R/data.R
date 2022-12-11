#' The data are derived from five real gut microbiome studies of colon cancer.
#' @format A list with 3 elements
#' \describe{
#'    \item{OTU}{a list contains OTU counts from five different studies}
#'    \item{Tax}{a matrix of taxonomy table from Rank1 (kingdom level) to Rank6 (genus level)}
#'    \item{covariate}{covariates of interest}
#' }
"data.meta"

#' Genus-level taxa count data for 574 observations.
#'
#' A dataset containing genus-level count data pooled from five studies
#'
#' @format A 574 x 133 Matrix
#' \describe{
#'    \item{rownames}{sample id}
#'    \item{colnames}{generic name for each taxa}
#' }
"count.genus"

#' Subject-level data for 574 observations.
#'
#' A dataset containing subject-level variables pooled from five gut microbiome studies of colon cancer.
#'
#' @format A dataframe
#' \describe{
#'    \item{Sample_id}{Sample id for each observation}
#'    \item{External_id}{External_id}
#'    \item{Age}{Age for each observation}
#'    \item{Gender}{Gender, F or M}
#'    \item{BMI}{BMI for each observation}
#'     ...
#'    \item{Group}{Sample belongs to case or control group}
#'     ...
#' }
"meta"

#' Taxonomy table
#'
#' Taxonomy ranks for the taxa counts in count.genus data
#'
#' @format A dataframe
#' \describe{
#'    \item{OTUID}{id}
#'    \item{kingdom}{the scientific name of the kingdom that the taxa belongs to}
#'    \item{phylum}{the scientific name of the phylum that the taxa belongs to}
#'    \item{class}{the scientific name of the class that the taxa belongs to}
#'    \item{order}{the scientific name of the order that the taxa belongs to}
#'    \item{family}{the scientific name of the family that the taxa belongs to}
#'    \item{genus}{generic name for each taxa conut}
#'    \item{mOTU}{mOTU}
#' }
"tax"
