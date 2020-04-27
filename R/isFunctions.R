#' ADAM classes checkers
#'
#' Functions to check if an object is of the specified class
#'
#' The list of functions includes:
#' \itemize{
#' \item \code{is.adam()} tests if the object was produced by a \link[adam]{adam} function
#' \item \code{is.adam.sim()} tests if the object was produced by \link[adam]{sim.adam} function;
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
is.adam <- function(x){
    return(inherits(x,"adam"))
}

#' @rdname isFunctions
#' @export
is.adam.sim <- function(x){
    return(inherits(x,"sim.adam"))
}
