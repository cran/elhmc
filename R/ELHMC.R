#'Empirical Likelihood Hamiltonian Monte Carlo Sampling
#'
#'This function draws samples from a Empirical Likelihood Bayesian posterior
#'distribution of parameters using Hamiltonian Monte Carlo.
#'
#'@param initial A vector containing the initial values of the parameter
#'@param data A matrix containing the data
#'@param fun The estimating function \eqn{g}. It takes in a parameter vector
#'\code{params} as the first argument and a data point vector \code{x} as the
#'  second parameter. This function returns a vector.
#'@param dfun A function that calculates the gradient of the estimating function
#'\eqn{g}. It takes in a parameter vector \code{params} as the first argument
#'  and a data point vector \code{x} as the second argument. This function
#'  returns a matrix.
#'
#'@param prior A function with one argument \code{x} that returns
#'  the prior density of the parameter of interest
#'@param dprior A function with one argument \code{x} that returns
#'  the gradient of the density of the parameter of interest
#'@param n.samples Number of samples to draw
#'@param lf.steps Number of leap frog steps in each Hamiltonian Monte Carlo
#'  update
#'@param epsilon The leap frog step size
#'@param tol EL tolerance
#'@param detailed If this is set to \code{TRUE}, the function will return a list
#'  with extra information.
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
#'  and the gradient of \eqn{g} as \code{dfun}.
#'@return The function returns a list with the following elements:
#'  \item{\code{samples}}{A matrix containing the parameter samples}
#'  \item{\code{acceptance.rate}}{The acceptance rate (numeric)}
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
#'## of the mean, using initial values (0.96, 0.97) and a standard normal distribution
#'## as the prior:
#'normal_prior <- function(x) {
#'  exp(-0.5 * x ^ 2) / sqrt(2 * pi)
#'}
#'normal_prior_gradient <- function(x) {
#'  -x * exp(-0.5 * x ^ 2) / sqrt(2 * pi)
#'}
#'set.seed(1234)
#'mean.samples <- ELHMC(initial = c(0.96, 0.97), data = x, fun = f, dfun = df,
#'                      n.samples = 50, prior = normal_prior,
#'                      dprior = normal_prior_gradient)
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
#'@export
#'
ELHMC <- function(initial, data, fun, dfun, prior, dprior,
                  n.samples = 100, lf.steps = 10, epsilon = 0.05, tol = 10^-5,
                  detailed = FALSE) {
  if(!(is.vector(initial) && is.numeric(initial))) {
    stop("initial must be a number or a numeric vector")
  }

  if(!(is.matrix(data) && is.numeric(data))) {
    stop("data must be a numeric matrix")
  }

  if(!CheckFuncArgs(fun, c("params", "x"))) {
    stop("fun must be a function with only two arguments \"params\" and \"x\"")
  }

  if(!CheckFuncArgs(dfun, c("params", "x"))) {
    stop("dfun must be a function with only two arguments \"params\" and \"x\"")
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

  if(epsilon <= 0) {
    stop("epsilon must be positive")
  }

  if(tol <= 0) {
    stop("tol must be positive")
  }

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
  for(i in 2:n.samples) {
    next.sample <- HMC(current.value, U = ELU, epsilon = epsilon,
                       lf.steps = lf.steps, detailed = detailed,
                       data = data, fun = fun, dfun = dfun,
                       prior = prior, dprior = dprior, tol = tol)
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
  }

  acceptance.rate <- sum(acceptance) / (n.samples - 1)

  if(!detailed) {
    return(list(samples = samples, acceptance.rate = acceptance.rate))
  }

  trajectory <- list(trajectory.q = trajectory.q,
                     trajectory.p = trajectory.p)
  results <- list(samples = samples,
                  acceptance.rate = acceptance.rate,
                  proposed = proposed,
                  acceptance = acceptance,
                  trajectory = trajectory)

  return(results)
}
