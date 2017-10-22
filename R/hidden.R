########## hidden functions to help in model implementation ##########

# parse formula and return a list that contains the model response
# matrix as element one, and the model matrix as element two
"parse.formula" <- function(formula, data = NULL) {
    # extract Y, X, and variable names for model formula and frame
    mf <- match.call(expand.dots = FALSE)
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf, "terms")

    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    xobs  <- dimnames(X)[[1]] # X observation names
    Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
    yname <- names(mf)[1]
    return(list(Y, yname, X, xvars))
}


# parse formula of bsam and return a list that contains the model response
# matrix and the model matrices of both parametric and nonparametric components
"interpret.bsam" <- function(formula) {
  # extract Y, X, and variable names for model formula and frame
  mf <- match.call(expand.dots = FALSE)
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  mt <- attr(mf, "terms")

  # null model support
  cvars <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
  cvars <- as.matrix(cvars) # matrix of predictors
  pnames <- dimnames(cvars)[[2]] # variable names of predictors
  xindex <- which(substr(pnames,1,3) == 'fs(')
  X <- cvars[, xindex, drop = FALSE]
  xnames <- gsub('[fs()]', '', colnames(X))
  colnames(X) <- xnames
  W <- cvars[, -xindex, drop = FALSE]
  wnames <- colnames(W)
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
  yname <- names(mf)[1]
  return(list(Y, yname, W, wnames, X, xnames))
}


# return a list that contains the shape constraints
# fmodel is a vector containing shape-types,
# fpm is a vector denoting increasing or decreasing
# nfun is the number of functions.
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
