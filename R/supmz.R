#' @name    supmz
#' @aliases supmz
#' @title Detecting Structural Change with Heteroskedasticity
#' @description Calculates the sup MZ value to detect the unknown structural break points under Heteroskedasticity
#'
#' @param formula Formula for the linear model to be used. It may contain any number of independent variables.
#' @param data    Data frame containing dependent and independent variables.
#' @param nBoot   Number of bootstrap samples to compute the critical region.
#'
#'
#' @return  MZ Gives values of MZ as given by Mumtaz et.al (2017)
#' @return  BreakLocation Provides the data point position where the structural break occured
#' @return  SupMzValue Returns the supremum value from MZ values
#' @return  SupMZ0 Returns the bootstrapped critical value for testing the significance of SupMZ
#' @return  nBoot Shows the number of bootstrap samples used to compute the critical region
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
#' @importFrom stats lm rnorm sigma quantile update
#' @import dplyr
#' @import magrittr
#'
#' @export
#'
#' @examples
#'
#' data(Japan)
#' fm1 <- supmz(formula = C~Y, data = Japan, nBoot = 10)
#' fm1
#'
#' data(Belgium)
#' fm2 <- supmz(formula = C~Y, data = Belgium, nBoot = 10)
#' fm2
#'
#' data(Srilanka)
#' fm3 <- supmz(formula = C~Y, data = Srilanka, nBoot = 10)
#' fm3
#'

supmz <- function(formula, data, nBoot = 100){
  UseMethod("supmz")
}
#' @export
#' @rdname supmz

supmz.default <- function(formula, data, nBoot = 100){
  k <- ifelse(test = dim(data)[1]<=10, yes = 3, no = 5)
  MZ     <- NULL
  fm0    <- lm(formula, data = data, subset = NULL)
  fitted <- as.numeric(fitted(fm0))
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
  n      <- which.max(MZ)
  SupMz  <- MZ[n]
  SupMZ01 <- NULL
  for(i in 1:nBoot){
    MZ0     <- NULL
    data1   <- data %>% mutate(ys=fitted + rnorm(n = nrow(data),mean = 0,sd = sqrt(sigma02)))
    fm00    <- update(fm0,ys~., data=data1)
    sigma00 <- sigma(fm00)^2
    for(V in k:(nrow(data1)-k)){
      subdata1 <- data1 %>% slice(1:V)
      fm10     <- update(fm1,ys~.,data= subdata1)
      sigma102 <- sigma(fm10)^2
      subdata2 <- data1 %>% slice(V+1:nrow(data1))
      fm20     <- update(fm2,ys~., data = subdata2)
      sigma202 <- sigma(fm20)^2
      MZ0[V+1] <- (nrow(data1)- length(fm10$coefficients))*log(sigma00,base = exp(1))-
        (V-length(fm10$coefficients))*log(sigma102,base = exp(1))-
        (nrow(data)-V-length(fm10$coefficients))*log(sigma202,base = exp(1))
    }
    n0        <- which.max(MZ0)
    SupMZ01[i] <- MZ0[n0]
  }

  return(
    list(
        MZ            = MZ
      , BreakLocation = n
      , SupMzValue    = SupMz
      , SupMZ0        = quantile(x = SupMZ01, probs = 0.95)
      , nBoot         = nBoot
      )
      )
}
