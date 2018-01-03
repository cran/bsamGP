intsim <- function(f, delta) {
  r <- length(f)
  if (r == 2 * floor(r/2)) {
    stop("ERROR: Even number of rows for simpson's integration")
  } else if (r == 3) {
    t <- c(1, 4, 1)
    fint <- sum(t * f) * delta/3
  } else {
    t <- c(1, rep(c(4, 2), times = (r - 3)/2), 4, 1)
    fint <- sum(t * f) * delta/3
  }
  fint
}
