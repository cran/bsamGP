"intgrat" <- function(f, delta){
  nr <- nrow(f)
  nc <- ncol(f)
  f1 <- matrix(t(f[1,]), nrow = nc, ncol = 1)
  fn <- matrix(t(f[nr,]), nrow = nc, ncol = 1)
  f <- matrix(f[2:(nr-1),], nrow = nr - 2, ncol = nc)
  fint <- delta*(apply(f, 2, sum) + (f1 + fn)/2)
  return(fint)
}