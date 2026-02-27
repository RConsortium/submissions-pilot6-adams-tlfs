# Load required packages
library(dplyr)        # Data manipulation

#' Derive pooled site group
#'
#' `format_sitegr1` derives the `SITEGR1` variable, which pools sites with fewer than 3 subjects in any one treatment
#' group.
#'
#' @param x The input vector of values to format
#'
#' @export
format_sitegr1 <- function(x) {
  case_when(
    x %in% c("702", "706", "707", "711", "714", "715", "717") ~ "900",
    TRUE ~ x
  )
}
