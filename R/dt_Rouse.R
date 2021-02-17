#'Rouse (1995) dataset
#'
#'Variables of the dataset is as follows:
#'\describe{
#' \item{educ86}{Years of education since 1986.}
#' \item{twoyr}{Attending a two-year college immediately after high school.}
#' \item{female}{Gender: 1 if female and 0 otherwise.}
#' \item{black}{Race: 1 if African American and 0 otherwise.}
#' \item{hispanic}{Race: 1 if Hispanic and 0 otherwise.}
#' \item{bytest}{Test score.}
#' \item{fincome}{Family income.}
#' \item{fincmiss}{Missingness indicator for family income.}
#' \item{IV}{Instrumental variable: encouagement to attend a two-year college.}
#' \item{dadeduc}{Dad's education: College - 2; Some college - 1; Neither - 0.}
#' \item{momeduc}{Mom's education: College - 2; Some college - 1; Neither - 0.}
#'}
#'
#'@docType data
#'
#'@usage data(dt_Rouse)
#'
#'@keywords datasets
#'
#'@format A data frame with 3037 rows, 8 observed variables, 1 binary instrumental variable,
#'  1 treatment, and 1 continuous response.
#'@source ss
"dt_Rouse"
