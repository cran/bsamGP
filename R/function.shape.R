"function.shape" <- function(shape = c("Free", "Increasing", "Decreasing", "IncreasingConvex", "DecreasingConcave", "IncreasingConcave",
                                       "DecreasingConvex", "IncreasingS", "DecreasingS", "IncreasingRotatedS", "DecreasingRotatedS", "InvertedU", "Ushape")) {
  choices <- c("Free", "Increasing", "Decreasing", "IncreasingConvex", "DecreasingConcave", "IncreasingConcave", "DecreasingConvex",
               "IncreasingS", "DecreasingS", "IncreasingRotatedS", "DecreasingRotatedS", "InvertedU", "Ushape")
  shape <- match.arg(shape, choices, several.ok = TRUE)
  nfun <- length(shape)
  fmodel <- numeric(nfun)
  fpm <- numeric(nfun)
  for (i in 1:nfun) {
    switch(shape[i], Free = {
      fmodel[i] <- 1
      fpm[i] <- 1
    }, Increasing = {
      fmodel[i] <- 2
      fpm[i] <- 1
    }, Decreasing = {
      fmodel[i] <- 2
      fpm[i] <- -1
    }, IncreasingConvex = {
      fmodel[i] <- 3
      fpm[i] <- 1
    }, DecreasingConcave = {
      fmodel[i] <- 3
      fpm[i] <- -1
    }, IncreasingConcave = {
      fmodel[i] <- 4
      fpm[i] <- 1
    }, DecreasingConvex = {
      fmodel[i] <- 4
      fpm[i] <- -1
    }, IncreasingS = {
      fmodel[i] <- 5
      fpm[i] <- 1
    }, DecreasingS = {
      fmodel[i] <- 5
      fpm[i] <- -1
    }, IncreasingRotatedS = {
      fmodel[i] <- 6
      fpm[i] <- 1
    }, DecreasingRotatedS = {
      fmodel[i] <- 6
      fpm[i] <- -1
    }, InvertedU = {
      fmodel[i] <- 7
      fpm[i] <- 1
    }, Ushape = {
      fmodel[i] <- 7
      fpm[i] <- -1
    })
  }
  list(fmodel = fmodel, fpm = fpm, nfun = nfun)
}
