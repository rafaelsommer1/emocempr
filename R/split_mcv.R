#' Takes a variable with multiple numeric variables and split into a data.frame with unique values in each col
#' @import dplyr
#' @import stringr
#' @importFrom utils read.csv
#' @author R.C.S
#' @param var Takes a string variable from a data frame with 'n' possible values
#' @param n Integer with the maximum number of possible values in a individual subject
#' @param dataf Data frame from which the variable will be taken
#' @return a data frame
#' @export



split_mcv <- function (var, n, dataf) {

  # garante que n é inteiro
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n <= 0) stop("Por favor insira um valor válido para n (inteiro > 0).")

  charc <- as.character(dataf[[var]])

  # Main
  split <- stringr::str_split_fixed(charc, ",", n)

  df <- as.data.frame(split, stringsAsFactors = FALSE)

  for (i in 1:n) {
    df[, i] <- as.numeric(df[, i])
    df[, i][is.na(df[, i])] <- 0
  }

  return(df)
}


