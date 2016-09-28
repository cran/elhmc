ELU <- function(x, data, fun, dfun, prior, dprior, tol) {
  fun.results <- plyr::aaply(data, 1, function(d) {
    fun(params = x, x = d)
  })

  fun.gradients <- plyr::aaply(data, 1, function(d) {
    dfun(params = x, x = d)
  })

  density <- prior(x)
  density.gradient <- dprior(x)

  el <- emplik::el.test(fun.results, rep(0, length(x)))

  if(abs(mean(el$wts) - 1) > tol) {
    CustomStop("invalid_el_weight_sum",
                "Sum of EL weights is not close enough to 1.")
  }

  u <- - sum(log(density)) - sum(log(el$wts / length(el$wts)))

  dellogL <- array(0, c(length(el$wts), length(x)))
  for (i in 1:NROW(data))
  {
    dellogL[i, ] <- el$wts[i] * t(as.matrix(el$lambda)) %*% fun.gradients[i, , ]
  }
  gradient <- apply(dellogL, 2, sum) - density.gradient

  result <- u
  attr(result, "gradient") <- gradient

  return(result)
}
