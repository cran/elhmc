#'@importFrom stats runif rnorm
HMC <- function(initial, U, epsilon = 0.05, lf.steps = 10, detailed = FALSE, ...) {
  current.q <- q <- initial
  current.p <- p <- rnorm(length(initial), 0, 1)

  if(detailed) {
    trajectory.q <- matrix(NA, lf.steps + 1, length(initial))
    trajectory.p <- matrix(NA, lf.steps + 2, length(initial))

    count.q <- 1
    count.p <- 1

    trajectory.q[count.q, ] <- q
    trajectory.p[count.p, ] <- p
  }

  result <- tryCatch({
    u <- current.U <- U(q, ...)
    p <- p - epsilon * attr(u, "gradient") / 2

    if(detailed) {
      count.p <- count.p + 1
      trajectory.p[count.p, ] <- p
    }

    if(lf.steps > 1) {
      for(i in 1:(lf.steps - 1)) {
        q <- q + epsilon * p

        if(detailed) {
          count.q <- count.q + 1
          trajectory.q[count.q, ] <- q
        }

        u <- U(q, ...)
        p <- p - epsilon * attr(u, "gradient")

        if(detailed) {
          count.p <- count.p + 1
          trajectory.p[count.p, ] <- p
        }
      }
    }

    q <- q + epsilon * p

    if(detailed) {
      count.q <- count.q + 1
      trajectory.q[count.q, ] <- q
    }

    u <- proposed.U <- U(q, ...)
    p <- - (p - epsilon * attr(u, "gradient") / 2)

    if(detailed) {
      count.p <- count.p + 1
      trajectory.p[count.p, ] <- p
    }

    current.K <- sum(current.p ^ 2) / 2
    proposed.K <- sum(p ^ 2) / 2

    proposed.value <- q
    accepted <- FALSE
    if(runif(1) < exp(current.U - proposed.U + current.K - proposed.K)) {
      accepted <- TRUE
    }

    if(detailed) {
      list(current.value = current.q,
           proposed.value = proposed.value,
           accepted = accepted,
           trajectory = list(trajectory.q = trajectory.q,
                             trajectory.p = trajectory.p))
    } else {
      list(current.value = current.q,
           proposed.value = proposed.value,
           accepted = accepted)
    }
  }, invalid_el_weight_sum = function(e) {
    if(detailed) {
      trajectory.q <- trajectory.q[1:count.q, ]
      trajectory.p <- trajectory.p[1:count.p, ]

      result <- list(current.value = current.q,
                     proposed.value = rep(NA, length(initial)),
                     accepted = FALSE,
                     trajectory = list(trajectory.q = trajectory.q,
                                       trajectory.p = trajectory.p))
    } else {
      result <- list(current.value = current.q,
                     proposed.value = rep(NA, length(initial)),
                     accepted = FALSE)
    }

    return(result)
  })

  return(result)
}
