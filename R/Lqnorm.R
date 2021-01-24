#' @title Lqnorm function
#' @description generate q's norm of a vector 
#' @param q integer 
#' @param v vector
#' @return a real number which is q's norm of a vector
#' @details q-norm (also called lq-norm)
#' @examples 
#' when q=2, it's called Euclidean norm
#' @rdname Lqnorm
#' @export
#' 
#'   
Lqnorm <- function(q, v){
  return((sum(abs(v) ^ q)) ^ (1 / q))
}