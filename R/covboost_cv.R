## Copyright (C) 2020 Berent Lunde
## License: GPL-3

#' Boosted Covariance Matrix Estimation With CV
#'
#' @param x An \code{n x p} numeric matrix
#' @param learning_rate Scaling the path of elements in the covariance matrix
#' @param niter The number of boosting iterations
#' @param nfolds The number of folds for k-fold cross validation
#' @param cores The number of cores use for parallel computations
#'
#' @return A list of useful cross validation information from boosting out the covariance matrix
#' @examples
#' \dontrun{
#' ## Generate data with random correlation
#' p <- 50
#' n <- 50
#' x <- matrix(nrow=n, ncol=p)
#' x[,1] <- rnorm(n)
#' for(j in 2:p){
#'     rho <- runif(1,-1,1)*rbinom(1,1,0.5)
#'     col_dependence <- sample(j-1,1)
#'     x[,j] <- rnorm(n, mean=rho*x[,col_dependence], sd=(1-rho^2))
#' }
#' lrn_rate <- 0.5
#' cov_cv <- covboost_cv(x, learning_rate=lrn_rate)
#' cov_cv$cvplot
#' cov_cv$opt_iter
#' }
#'
#' @importFrom stats cov cov2cor quantile
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @rdname covboost_cv
#' @export
covboost_cv <- function(x, learning_rate=0.5, niter=1000, nfolds=10, cores=1)
{
    # Boosts out a covariance matrix from the Identity matrix
    # for p>>n matrix will become singular. The function notices this and terminates
    # Optimal covariance matrix is likely reached some iterations before this

    # requires: ggplot2

    # x: data
    # learning_rate: shrink each step
    # niter: max number of iterations
    # nfolds: number of folds in CV
    # cores: number of cores used in parallel

    `10-fold cv` <- cv <- Var1 <- Var2 <- iterations <- qlower <- qupper <- value <- NULL

    # Check input
    error_messages <- c()
    error_messages_type <- c(
        "x" = "\n Error: x must be a numeric matrix",
        "lrn_rate" = "\n Error: learning_rate must be a number between 0 and 1",
        "niter" = "\n Error: niter must be an integer >= 1",
        "nfolds" = "\n Error: nfolds must be an integer > 1 and <= nrow(x)",
        "cores" = "\n Error: cores must be an integer >= 1"
    )

    # check x
    if(!is.matrix(x) || !is.numeric(x)){
        error_messages <- c(error_messages, error_messages_type["x"])
    }

    # learning_rate
    if(!is.numeric(learning_rate) ||
       !(length(learning_rate)==1) ||
       learning_rate < 0 ||
       learning_rate > 1){
        error_messages <- c(error_messages, error_messages_type["lrn_rate"])
    }

    # niter
    if(!is.numeric(niter) ||
       length(niter)>1 ||
       niter < 1){
        error_messages <- c(error_messages, error_messages_type["niter"])
    }

    # nfolds
    if(!is.numeric(nfolds) ||
       length(nfolds)>1 ||
       nfolds < 1 ||
       nfolds > nrow(x)){
        error_messages <- c(error_messages, error_messages_type["nfolds"])
    }

    # cores
    if(!is.numeric(cores) ||
       length(cores)>1 ||
       cores < 1){
        error_messages <- c(error_messages, error_messages_type["cores"])
    }

    # Any error messages?
    if(length(error_messages)>0)
        stop(error_messages)


    # Dimensions
    n <- nrow(x)
    p <- ncol(x)

    # Set up variables
    MAX_CONSECUTIVE_NO_IMPROVEMENT <- 50
    NO_IMPROVEMENT_COUNTER <- 0
    BEST <- NULL
    cv_nll_cor_k <- numeric(niter+1)
    cv_nll_cor <- matrix(nrow=niter+1, ncol=nfolds)
    stops <- rep(niter+1, nfolds)

    pb <- txtProgressBar(min=0, max=(niter+1)*nfolds, style=3)
    for(k in 1:nfolds)
    {
        # reset counter
        NO_IMPROVEMENT_COUNTER <- 0

        # define holdout set
        holdout <- split(sample(n,n), 1:nfolds)

        # starting matrices
        Ak <- cov(x[-holdout[[k]],])
        Vk <- diag(diag(Ak))
        eDiag <- sqrt(Vk)
        invDiag <- sqrt(diag(1/diag(Vk)))
        Ak_rho <- invDiag %*% Ak %*% invDiag
        Bk_rho <- diag(p) # Identity: inverse of diagonal variance matrix
        Bk <- round(eDiag %*% Bk_rho %*% eDiag, digits=9) # = Vk = diag(diag(Ak))
        cv_nll_cor_k[1] <- NA
        cv_nll_cor_k[1] <- try(-sum(.dmvnorm_arma_mc(x[holdout[[k]],], rep(0,p), Bk, logd = TRUE, cores=cores)), silent = T)
        BEST <- cv_nll_cor_k[1]

        for(i in 2:(niter+1))
        {
            # difference
            Dk_rho <- Ak_rho - Bk_rho
            # COOL NOTE: POSSIBLE TO HELP WITH GRAPH!!
            ind <- which(abs(Dk_rho)==max(abs(Dk_rho)), arr.ind=T)

            # update Bk_rho
            Bk_rho[ind] <- Bk_rho[ind] + learning_rate*Dk_rho[ind]

            # check holdout loss
            Bk <- round(eDiag %*% Bk_rho %*% eDiag, digits=9) # transform to correlation, round ensures symmetry
            cvnll_i_k <- NA # default if fails
            cvnll_i_k <- try(-sum(.dmvnorm_arma_mc(x[holdout[[k]],], rep(0,p), Bk, logd = TRUE, cores=cores)), silent = T)

            # checks
            if(!is.finite(cvnll_i_k)) {
                # fill remaining
                cv_nll_cor_k[i:(niter+1)] <- cv_nll_cor_k[i-1]
                # get stop-point
                stops[k] <- i-1
                #cat("stopping at iteration: ", i, "\n",
                #    "updating niter to ", i-1)
                break
            }else{
                #update
                cv_nll_cor_k[i] <- cvnll_i_k

                # check counter
                if(cv_nll_cor_k[i]>=BEST){
                    NO_IMPROVEMENT_COUNTER <- NO_IMPROVEMENT_COUNTER + 1
                }else{
                    NO_IMPROVEMENT_COUNTER <- 0
                }
            }

            # check counter
            if(NO_IMPROVEMENT_COUNTER==MAX_CONSECUTIVE_NO_IMPROVEMENT){
                cv_nll_cor_k[i:(niter+1)] <- cv_nll_cor_k[i-1]
                stops[k] <- i
                break
            }

            setTxtProgressBar(pb, value=(k-1)*niter+i)
        }

        cv_nll_cor[,k] <- cv_nll_cor_k

    }
    close(pb)

    cat("preparing data...\n")
    cv_nll <- cv_nll_cor[1:max(stops),,drop=FALSE]
    cv_nll_mean <- rowMeans(cv_nll, na.rm=T)
    cv_nll_qupper <- sapply(1:nrow(cv_nll), function(i){quantile(cv_nll[i,], 0.6, na.rm = T)})
    cv_nll_qlower <- sapply(1:nrow(cv_nll), function(i){quantile(cv_nll[i,], 0.4, na.rm=T)})
    chol_decomp_fails <- stops[which(stops < (niter+1))]

    # update niter and create data
    niter <- nrow(cv_nll)
    .df <- data.frame(1:niter, cv_nll_mean[1:niter])
    .df2 <- data.frame(1:niter, cv_nll_qupper[1:niter], cv_nll_qlower[1:niter])
    names(.df) <- c("iterations", "10-fold cv")
    names(.df2) <- c("iterations", "qupper", "qlower")

    # CV loss vs iterations plot
    cat("preparing loss...\n")
    cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbp2) <- colnames(.df)
    opt_iter = .df[which.min(.df$`10-fold cv`), "iterations"]

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        .p <- ggplot2::ggplot(data=.df) +
            ggplot2::geom_point(ggplot2::aes(x=iterations, y=`10-fold cv`), colour=cbp2[2]) +
            ggplot2::geom_line(ggplot2::aes(x=iterations, y=`10-fold cv`), colour=cbp2[2]) +
            ggplot2::geom_ribbon(data=.df2, ggplot2::aes(iterations, ymax=qupper, ymin=qlower), fill=cbp2[2], alpha=0.5) +
            ggplot2::geom_vline(xintercept=opt_iter, size=1) +
            ggplot2::xlab("Iteration") +
            ggplot2::ylab("Gaussian loss") +
            ggplot2::ggtitle("Gaussian CV Loss Versus Iterations") +
            ggplot2::theme_bw() +
            ggplot2::geom_vline(xintercept=chol_decomp_fails, size=1, colour="blue")
    } else {
        .p <- cat("Install ggplo2 package to get cv-plots \n")
    }

    cat("preparing results...\n")
    .res_data <- data.frame(1:niter, cv_nll_mean[1:niter], cv_nll_qupper[1:niter], cv_nll_qlower[1:niter])
    names(.res_data) <- c("iterations", "cv-mean", "cv_6_quantile", "cv_4_quantile")
    res <- list(
        #opt_shrinkage = d_min_loss,
        opt_iter = opt_iter,
        cvplot = .p,
        data = .res_data,
        chol_decomp_fails = chol_decomp_fails
    )


    return(res)

}
