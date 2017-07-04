#'Empirical Likelihood Hamiltonian Monte Carlo Sampling
#'
#'This function draws samples from a Empirical Likelihood Bayesian posterior
#'distribution of parameters using Hamiltonian Monte Carlo.
#'
#'@param initial a vector containing the initial values of the parameters
#'@param data a matrix containing the data
#'@param fun the estimating function \eqn{g}. It takes in a parameter vector
#'\code{params} as the first argument and a data point vector \code{x} as the
#'  second parameter. This function returns a vector.
#'@param dfun a function that calculates the gradient of the estimating function
#'\eqn{g}. It takes in a parameter vector \code{params} as the first argument
#'  and a data point vector \code{x} as the second argument. This function
#'  returns a matrix.
#'@param prior a function with one argument \code{x} that returns the prior
#'  densities of the parameters of interest
#'@param dprior a function with one argument \code{x} that returns
#'  the gradients of the log densities of the parameters of interest
#'@param n.samples number of samples to draw
#'@param lf.steps number of leap frog steps in each Hamiltonian Monte Carlo
#'  update
#'@param epsilon the leap frog step size(s). This has to be a single numeric
#'  value or a vector of the same length as \code{initial}.
#'@param p.variance the covariance matrix of a multivariate normal distribution
#'  used to generate the initial values of momentum \eqn{p} in
#'  Hamiltonian Monte Carlo. This can also be a single numeric value or
#'  a vector. See Details.
#'@param tol EL tolerance
#'@param detailed If this is set to \code{TRUE}, the function will return a list
#'  with extra information.
#'@param FUN the same as \code{fun} but takes in a matrix \code{X} instead of
#'  a vector \code{x} and returns a matrix so that \code{FUN(params, X)[i, ]} is
#'  the same as \code{fun(params, X[i, ])}. Only one of \code{FUN} and
#'  \code{fun} should be provided. If both are then \code{fun} is ignored.
#'@param DFUN the same as \code{dfun} but takes in a matrix \code{X} instead of
#'  a vector \code{x} and returns an array so that \code{DFUN(params, X)[, , i]}
#'  is the same as \code{dfun(params, X[i, ])}. Only one of \code{DFUN} and
#'  \code{dfun} should be provided. If both are then \code{dfun} is ignored.
#'@details Suppose there are data \eqn{x = (x_1, x_2, ..., x_n)} where \eqn{x_i}
#'  takes values in \eqn{R^p} and follow probability distribution \eqn{F}.
#'  Also, \eqn{F} comes from a family of distributions that depends on
#'  a parameter \eqn{\theta = (\theta_1, ..., \theta_d)} and there is
#'  a smooth function
#'  \eqn{g(x_i, \theta) = (g_1(x_i, \theta), ...,g_q(x_i, \theta))^T} that
#'  satisfies \eqn{E_F[g(x_i, \theta)] = 0} for \eqn{i = 1, ...,n}.
#'
#'  \code{ELHMC} draws samples from a Empirical Likelihood Bayesian
#'  posterior distribution of the parameter \eqn{\theta}, given the data \eqn{x}
#'  as \code{data}, the smoothing function \eqn{g} as \code{fun},
#'  and the gradient of \eqn{g} as \code{dfun} or \eqn{G(X) =
#'  (g(x_1), g(x_2), ..., g(x_n))^T} as \code{FUN} and the gradient of \eqn{G}
#'  as \code{DFUN}.
#'@return The function returns a list with the following elements:
#'  \item{\code{samples}}{A matrix containing the parameter samples}
#'  \item{\code{acceptance.rate}}{The acceptance rate}
#'  \item{\code{call}}{The matched call}
#'
#'  If \code{detailed = TRUE}, the list contains these extra elements:
#'  \item{\code{proposed}}{A matrix containing the proposed values at
#'    \code{n.samaples - 1} Hamiltonian Monte Carlo updates}
#'  \item{\code{acceptance}}{A vector of \code{TRUE}/\code{FALSE} values
#'    indicates whether each proposed value is accepted}
#'  \item{\code{trajectory}}{A list with 2 elements \code{trajectory.q} and
#'    \code{trajectory.p}. These are lists of matrices contraining position and
#'    momentum values along trajectory in each Hamiltonian Monte Carlo update.}
#'@examples
#'\dontrun{
#'## Suppose there are four data points (1, 1), (1, -1), (-1, -1), (-1, 1)
#'x = rbind(c(1, 1), c(1, -1), c(-1, -1), c(-1, 1))
#'## If the parameter of interest is the mean, the smoothing function and
#'## its gradient would be
#'f <- function(params, x) {
#'  x - params
#'}
#'df <- function(params, x) {
#'  rbind(c(-1, 0), c(0, -1))
#'}
#'## Draw 50 samples from the Empirical Likelihood Bayesian posterior distribution
#'## of the mean, using initial values (0.96, 0.97) and standard normal distributions
#'## as priors:
#'normal_prior <- function(x) {
#'    exp(-0.5 * x[1] ^ 2) / sqrt(2 * pi) * exp(-0.5 * x[2] ^ 2) / sqrt(2 * pi)
#'}
#'normal_prior_log_gradient <- function(x) {
#'    -x
#'}
#'set.seed(1234)
#'mean.samples <- ELHMC(initial = c(0.96, 0.97), data = x, fun = f, dfun = df,
#'                      n.samples = 50, prior = normal_prior,
#'                      dprior = normal_prior_log_gradient)
#'plot(mean.samples$samples, type = "l", xlab = "", ylab = "")
#'}
#'@references Chaudhuri, S., Mondal, D. and Yin, T. (2015)
#'  Hamiltonian Monte Carlo sampling in Bayesian empirical likelihood
#'  computation.
#'  \emph{Journal of the Royal Statistical Society: Series B}.
#'
#'  Neal, R. (2011) MCMC for using Hamiltonian dynamics.
#'  \emph{Handbook of Markov Chain Monte Carlo}
#'  (eds S. Brooks, A.Gelman, G. L.Jones and X.-L. Meng), pp. 113-162.
#'  New York: Taylor and Francis.
#'@importFrom plyr aaply
#'@export
#'
ELHMC <- function(initial, data, fun, dfun, prior, dprior,
                  n.samples = 100, lf.steps = 10, epsilon = 0.05,
                  p.variance = 1, tol = 10^-5,
                  detailed = FALSE, FUN, DFUN) {
  if(!(is.vector(initial) && is.numeric(initial))) {
    stop("initial must be a number or a numeric vector")
  }

  if(!(is.matrix(data) && is.numeric(data))) {
    stop("data must be a numeric matrix")
  }
  
  if(!CheckFuncArgs(prior, "x")) {
    stop("prior must be a function with only one argument \"x\"")
  }

  if(!CheckFuncArgs(dprior, "x")) {
    stop("dprior must be a function with only one argument \"x\"")
  }

  if(n.samples <= 1) {
    stop("n.samples must be larger than 1")
  }

  if(lf.steps < 1) {
    stop("lf.steps must be at least 1")
  }
  
  if(!(is.vector(epsilon) && is.numeric(epsilon) &&
       (length(epsilon) == 1) || length(epsilon) == length(initial))) {
    stop(paste("epsilon must be a single numeric value or",
               "a numeric vector of the same length as initial"))
  }

  if(any(epsilon <= 0)) {
    stop("epsilon must be all positive")
  }
  
  if(!((is.vector(p.variance) && is.numeric(p.variance) &&
       (length(p.variance) == 1) || length(p.variance) == length(initial)) ||
       (is.matrix(p.variance) && ncol(p.variance) == nrow(p.variance) &&
        ncol(p.variance) == length(initial)))) {
    stop(paste("p.variance must be a single numeric value or",
               "a numeric vector of the same length as initial or",
               "a numeric square matrix with number of rows equal to",
               "length(initial)"))
  }
  
  if(any(p.variance < 0)) {
    stop("p.variance must be all at least 0")
  }

  if(tol <= 0) {
    stop("tol must be positive")
  }
  
  if(missing(FUN)) {
    if(missing(fun)) {
      stop("either fun or FUN (but not both) must be provided")
    }
    
    if(!CheckFuncArgs(fun, c("params", "x"))) {
      stop("fun must be a function with only two arguments \"params\" and \"x\"")
    }
    
    FUN <- function(params, X) {
      unname(aaply(X, 1, function(d) {
        fun(params = params, x = d)
      }, .drop = FALSE))
    }
  } else {
    if(!CheckFuncArgs(FUN, c("params", "X"))) {
      stop("FUN must be a function with only two arguments \"params\" and \"X\"")
    }
  }
  
  if(missing(DFUN)) {
    if(missing(dfun)) {
      stop("either dfun or DFUN (but not both) must be provided")
    }
    
    if(!CheckFuncArgs(dfun, c("params", "x"))) {
      stop("dfun must be a function with only two arguments \"params\" and \"x\"")
    }
    
    DFUN <- function(params, X) {
      result <- unname(aaply(X, 1, function(d) {
        dfun(params = params, x = d)
      }, .drop = FALSE))
      result <- unname(aperm(result, c(2, 3, 1)))
    }
  } else {
    if(!CheckFuncArgs(DFUN, c("params", "X"))) {
      stop("DFUN must be a function with only two arguments \"params\" and \"X\"")
    }
  }
  
  cl <- match.call()

  n.samples = floor(n.samples)
  lf.steps = floor(lf.steps)

  samples <- matrix(NA, nrow = n.samples, ncol = length(initial))
  samples[1, ] <- initial
  acceptance <- rep(NA, n.samples - 1)

  if(detailed) {
    proposed <- matrix(NA, nrow = n.samples - 1, ncol = length(initial))

    trajectory.q <- list()
    length(trajectory.q) <- n.samples - 1

    trajectory.p <- list()
    length(trajectory.p) <- n.samples - 1
  }

  current.value <- initial
  
  progress.bar <- utils::txtProgressBar(min = 1, max = n.samples, initial = 1,
                                        style = 3)
  on.exit(close(progress.bar))
  
  for(i in 2:n.samples) {
    next.sample <- HMC(current.value, U = ELU, epsilon = epsilon,
                       lf.steps = lf.steps, detailed = detailed,
                       data = data, fun = FUN, dfun = DFUN,
                       prior = prior, dprior = dprior, p.variance = p.variance,
                       tol = tol)
    samples[i, ] <- current.value <- if(next.sample$accepted) {
      next.sample$proposed.value
    } else {
      next.sample$current.value
    }
    acceptance[i - 1] <- next.sample$accepted

    if(detailed) {
      proposed[i - 1, ] <- next.sample$proposed.value

      trajectory.q[[i - 1]] <- next.sample$trajectory$trajectory.q
      trajectory.p[[i - 1]] <- next.sample$trajectory$trajectory.p
    }
    
    utils::setTxtProgressBar(progress.bar, i)
  }

  acceptance.rate <- sum(acceptance) / (n.samples - 1)

  if(!detailed) {
    return(list(samples = samples, acceptance.rate = acceptance.rate,
                call = cl))
  }

  trajectory <- list(trajectory.q = trajectory.q,
                     trajectory.p = trajectory.p)
  results <- list(samples = samples,
                  acceptance.rate = acceptance.rate,
                  proposed = proposed,
                  acceptance = acceptance,
                  trajectory = trajectory,
                  call = cl)

  return(results)
}
