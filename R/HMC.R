#'@importFrom stats runif rnorm
#'@importFrom MASS mvrnorm
HMC <- function(initial, U, epsilon = 0.05, lf.steps = 10, p.variance = 1,
                detailed = FALSE, ...) {
  n <- length(initial)
  
  if(!is.matrix(p.variance)) {
    if(length(p.variance) == 1) {
      p.variance <- diag(rep(p.variance, n))
    } else {
      p.variance <- diag(p.variance)
    }
  }
  
  current.q <- q <- initial
  current.p <- p <- mvrnorm(1, rep(0, n), p.variance)

  if(detailed) {
    n.row.q <- lf.steps + 1
    n.row.p <- lf.steps + 2
    
    trajectory.q <- matrix(NA, n.row.q, n)
    trajectory.p <- matrix(NA, n.row.p, n)

    count.q <- 1
    count.p <- 1

    trajectory.q[count.q, ] <- q
    trajectory.p[count.p, ] <- p
    
    step.q <- 0
    step.p <- 0
    
    rownames(trajectory.q) <- rep(0, n.row.q)
    rownames(trajectory.p) <- rep(0, n.row.p)
  }

  result <- tryCatch({
    u <- current.U <- U(q, ...)
    p <- p - epsilon * attr(u, "gradient") / 2

    if(detailed) {
      count.p <- count.p + 1
      trajectory.p[count.p, ] <- p
      step.p <- step.p + 0.5
      rownames(trajectory.p)[count.p] <- step.p
    }

    if(lf.steps > 1) {
      for(i in 1:(lf.steps - 1)) {
        q <- q + epsilon * p

        if(detailed) {
          count.q <- count.q + 1
          trajectory.q[count.q, ] <- q
          step.q <- step.q + 1
          rownames(trajectory.q)[count.q] <- step.q
        }

        u <- U(q, ...)
        p <- p - epsilon * attr(u, "gradient")

        if(detailed) {
          count.p <- count.p + 1
          trajectory.p[count.p, ] <- p
          step.p <- step.p + 1
          rownames(trajectory.p)[count.p] <- step.p
        }
      }
    }

    q <- q + epsilon * p

    if(detailed) {
      count.q <- count.q + 1
      trajectory.q[count.q, ] <- q
      step.q <- step.q + 1
      rownames(trajectory.q)[count.q] <- step.q
    }

    u <- proposed.U <- U(q, ...)
    p <- - (p - epsilon * attr(u, "gradient") / 2)

    if(detailed) {
      count.p <- count.p + 1
      trajectory.p[count.p, ] <- p
      step.p <- step.p + 0.5
      rownames(trajectory.p)[count.p] <- step.p
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
  }, invalid_el_weight_mean = function(e) {
    if(detailed) {
      trajectory.q <- trajectory.q[1:count.q, ]
      trajectory.p <- trajectory.p[1:count.p, ]

      result <- list(current.value = current.q,
                     proposed.value = rep(NA, n),
                     accepted = FALSE,
                     trajectory = list(trajectory.q = trajectory.q,
                                       trajectory.p = trajectory.p))
    } else {
      result <- list(current.value = current.q,
                     proposed.value = rep(NA, n),
                     accepted = FALSE)
    }

    return(result)
  })

  return(result)
}
