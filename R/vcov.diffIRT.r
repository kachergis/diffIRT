vcov.diffIRT = function (object, ...) {
    if (is.null(object$hessian))
        stop("to obtain the covariance matrix of the parameter estimates, you should re-fit the model using 'se = TRUE'.\n")
    # the stored hessian is of -2LL; the covariance matrix is the inverse of the hessian of -LL
    covmat <- solve(.5 * object$hessian)
    covmat
}
