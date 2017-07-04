condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call),
    ...
  )
}

CustomStop <- function(subclass, message, call = sys.call(-1),
                        ...) {
  c <- condition(c(subclass, "error"), message, call = call, ...)
  stop(c)
}

CheckFuncArgs <- function(func, arguments) {
  if(!inherits(func, 'function')) {
    return(FALSE)
  }

  func.args <- formals(func)

  if(length(func.args) != length(arguments)) {
    return(FALSE)
  }

  if(!all(names(func.args) == arguments)) {
    return(FALSE)
  }

  return(TRUE)
}
