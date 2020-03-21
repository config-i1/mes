#' MES classes checkers
#'
#' Functions to check if an object is of the specified class
#'
#' The list of functions includes:
#' \itemize{
#' \item \code{is.mes()} tests if the object was produced by a \link[mes]{mes} function
#' \item \code{is.mes.sim()} tests if the object was produced by \link[mes]{sim.mes} function;
#' }
#'
#' @param x The object to check.
#' @return \code{TRUE} if this is the specified class and \code{FALSE} otherwise.
#'
#' @template ssAuthor
#' @keywords ts univar
#' @examples
#'
#' \dontrun{ourModel <- msarima(rnorm(100,100,10))}
#'
#' @rdname isFunctions
#' @export
is.mes <- function(x){
    return(inherits(x,"mes"))
}

#' @rdname isFunctions
#' @export
is.mes.sim <- function(x){
    return(inherits(x,"sim.mes"))
}
