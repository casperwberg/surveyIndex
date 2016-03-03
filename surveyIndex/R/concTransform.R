##' Concentration transform
##'
##' @title Helper function for plotting survey indices.
##' @param x a vector of log-responses
##' @return vector of transformed responses
concTransform <-
function (x)
{
    i <- order(x)
    ys <- sort(exp(x))
    p <- ys/sum(ys)
    x[i] <- cumsum(p)
    x
}
