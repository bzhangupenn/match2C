#' RHC dataset (after cleaning)
#'
#' Variables in the dataset are as follows.
#'
#' @format A data frame with 2,998 rows and 27 variables:
#' \describe{
#'   \item{age_group}{Categorical age group (integer: 1, 2, 3): 1 = 18–29, 2 = 30–50, 3 = 51–65.}
#'   \item{edu}{Years of education.}
#'   \item{sexMale}{Binary indicator: 1 = male, 0 = female.}
#'   \item{raceother}{Binary indicator: 1 = race classified as Other, 0 otherwise.}
#'   \item{racewhite}{Binary indicator: 1 = White, 0 otherwise.  (If \code{raceother = 0} and \code{racewhite = 0}, race is Black.)}
#'   \item{income$25-$50k}{Binary indicator: household income between \$25k and \$50k.}
#'   \item{income> $50k}{Binary indicator: household income greater than \$50k.}
#'   \item{incomeUnder $11k}{Binary indicator: household income less than \$11k.
#'         (If \code{income$25-$50k = 0}, \code{income> $50k = 0}, and \code{incomeUnder $11k = 0}, income is \$11k–\$25k.)}
#'   \item{das2d3pc}{Duke Activity Status Index (DASI) score.}
#'   \item{caNo}{Binary indicator: 1 = no cancer, 0 otherwise.}
#'   \item{caYes}{Binary indicator: 1 = cancer present, 0 otherwise.  (If \code{caNo = 0} and \code{caYes = 0}, metastatic cancer.)}
#'   \item{temp1}{Body temperature.}
#'   \item{resp1}{Respiratory rate.}
#'   \item{paco21}{Arterial partial pressure of carbon dioxide, PaCO\eqn{_2}.}
#'   \item{ninsclasMedicare}{Binary indicator: Medicare only.}
#'   \item{ninsclasMedicare & Medicaid}{Binary indicator: dual Medicare/Medicaid.}
#'   \item{ninsclasNo insurance}{Binary indicator: no health insurance.}
#'   \item{ninsclasPrivate}{Binary indicator: private insurance only.}
#'   \item{ninsclasPrivate & Medicare}{Binary indicator: dual Private/Medicare.
#'         (If \code{ninsclasMedicare = 0}, \code{ninsclasMedicare & Medicaid = 0},
#'          \code{ninsclasNo insurance = 0}, \code{ninsclasPrivate = 0}, and
#'          \code{ninsclasPrivate & Medicare = 0}, insurance is Medicaid only.)}
#'   \item{wblc1}{White blood cell count.}
#'   \item{sod1}{Sodium.}
#'   \item{pot1}{Potassium.}
#'   \item{urin1adjusted}{Adjusted urine output.}
#'   \item{urin1miss}{Binary indicator: 1 if urine output data are missing, 0 otherwise.}
#'   \item{renalhx}{Binary indicator: 1 if the patient has a history of renal disease, 0 otherwise.}
#'   \item{liverhx}{Binary indicator: 1 if the patient has a history of liver disease, 0 otherwise.}
#'   \item{z}{Treatment indicator.}
#' }
#'
#' @docType data
#' @usage data(dt_rhc_ac)
#' @keywords datasets
#' @source See variable definitions at \url{https://hbiostat.org/data/repo/rhc}
"dt_rhc_ac"

#' RHC dataset (raw)
#'
#' Variables in the dataset are as follows.
#'
#' @format A data frame with 2,998 rows and 78 variables:
#' @docType data
#' @usage data(dt_rhc)
#' @keywords datasets
#' @source See variable definitions at \url{https://hbiostat.org/data/repo/rhc}
"dt_rhc"
