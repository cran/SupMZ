#' @name    supmz
#' @aliases supmz
#' @title Detecting Structural Change with Heteroskedasticity
#' @description Calculates the sup MZ value to detect the unknown structural break points under Heteroskedasticity
#'
#' @param formula Formula for the linear model to be used. It may contain any number of independent variables.
#' @param data    Data frame containing dependent and independent variables.
#'
#'
#' @return  MZ values as per given by Mumtaz et.al (2017)
#' @return  SupMz.Value i.e. the supremum value from MZ values.
#' @return  Break.Point.Location i.e. the data point number where the structural break occured.
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Sami Ullah      (\email{samiullahuos@@gmail.com})
#'          \item Gulfam Haider   (\email{haider.gulfam786@@gmail.com})
#'          }
#'
#' @references
#' Mumtaz Ahmed, Gulfam Haider & Asad Zaman (2017).
#' Detecting structural change with heteroskedasticity.
#' \emph{Communications in Statistics - Theory and Methods}.
#' \strong{46}(21):10446-10455,
#' DOI: 10.1080/03610926.2016.1235200
#'
#' @importFrom stats lm sigma
#' @import dplyr
#' @import magrittr
#'
#' @export
#'
#' @examples
#'
#' data(Japan)
#' fm1 <- supmz(formula = C~Y, data = Japan)
#' fm1
#'
#' data(Belgium)
#' fm2 <- supmz(formula = C~Y, data = Belgium)
#' fm2
#'
#' data(Srilanka)
#' fm3 <- supmz(formula = C~Y, data = Srilanka)
#' fm3
#'

supmz <- function(formula, data){
  UseMethod("supmz")
}
#' @export
#' @rdname supmz

supmz.default <- function(formula, data){
  k <- ifelse(test = dim(data)[1]<=10, yes = 3, no = 5)
  MZ  <- NULL
  fm0     <- lm(formula, data = data, subset = NULL)
  sigma02 <- sigma(fm0)^2
  for(V in k:(nrow(data)-k)){
    subdata1 <- data %>% slice(1:V)
    fm1      <- lm(formula, data = subdata1)
    sigma12  <- sigma(fm1)^2
    subdata2 <- data %>% slice(V+1:nrow(data))
    fm2      <- lm(formula, data = subdata2)
    sigma22  <- sigma(fm2)^2
    MZ[V+1]  <- (nrow(data)- length(fm0$coefficients))*log(sigma02,base = exp(1))-(V-length(fm0$coefficients))*log(sigma12,base = exp(1))-(nrow(data)-V-length(fm0$coefficients))*log(sigma22,base = exp(1))
    }
  n     <- which.max(MZ)
  SupMz <- MZ[n]
  return(list(MZ = MZ, SupMz.Value = SupMz, Break.Point.Location = n ))
}
