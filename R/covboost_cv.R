## Copyright (C) 2020 Berent Lunde
## License: GPL-3

#' Add together two numbers
#'
#' @param x An \code{n x p} matrix or data frame of dimensions n
#' @param learning_rate Scaling the path of elements in the covariance matrix
#' @param niter The number of boosting iterations
#' @param nfolds The number of folds for k-fold cross validation
#'
#' @return A list of useful cross validation information from boosting out the covariance matrix
#' @examples
#' add(1, 1)
#' add(10, 1)
covboost_cv <- function(x, learning_rate=0.01, niter=1000, nfolds=10)
{
    # Boosts out a covariance matrix from the Identity matrix
    # for p>>n matrix will become singular. The function notices this and terminates
    # Optimal covariance matrix is likely reached some iterations before this

    # requires: mvtnorm, ggplot2

    # x: data
    # learning_rate: shrink each step
    # niter: max number of iterations
    # nfolds: number of folds in CV

    n <- nrow(x)
    p <- ncol(x)

    cv_nll_k <- numeric(niter)
    cv_nll <- matrix(nrow=niter, ncol=nfolds)
    stops <- rep(niter, nfolds)

    cat("starting boosting...\n")
    pb <- txtProgressBar(min=0, max=niter*k, style=3)
    for(k in 1:nfolds)
    {
        Bk <- diag(p)
        holdout <- split(sample(n,n), 1:nfolds)

        for(i in 1:niter)
        {
            Ak <- cov(x[-holdout[[k]],])

            Dk <- Ak - Bk
            ind <- which(abs(Dk)==max(abs(Dk)), arr.ind=T)

            Bk[ind] <- Bk[ind] + learning_rate*Dk[ind]

            cvnll_i_k <- -sum(dmvnorm(x[holdout[[k]],], rep(0,p), Bk, log = TRUE))

            # checks
            if(!is.finite(cvnll_i_k)) {
                stops[k] <- i-1
                #cat("stopping at iteration: ", i, "\n",
                #    "updating niter to ", i-1)
                break
            }else{
                #update
                cv_nll_k[i] <- cvnll_i_k
            }

            setTxtProgressBar(pb, value=(k-1)*niter+i)

        }

        cv_nll[,k] <- cv_nll_k


    }
    close(pb)

    cat("preparing data...\n")
    # cut of at min stops
    cv_nll <- cv_nll[1:max(stops),]
    cv_nll_mean <- rowMeans(cv_nll, na.rm=T)
    cv_nll_qupper <- sapply(1:nrow(cv_nll), function(i){quantile(cv_nll[i,], 0.6, na.rm = T)})
    cv_nll_qlower <- sapply(1:nrow(cv_nll), function(i){quantile(cv_nll[i,], 0.4, na.rm=T)})

    # update niter and create data
    niter <- nrow(cv_nll)
    df <- data.frame(1:niter, cv_nll_mean[1:niter])
    df2 <- data.frame(1:niter, cv_nll_qupper[1:niter], cv_nll_qlower[1:niter])
    names(df) <- c("iterations", "10-fold cv")
    names(df2) <- c("iterations", "qupper", "qlower")

    # CV loss vs iterations plot
    cat("preparing loss...\n")
    cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    names(cbp2) <- colnames(df)
    opt_iter = df[which.min(df$`10-fold cv`), "iterations"]
    p <- df %>% #gather(type, `gaussian nll`, `10-fold cv`, factor_key = TRUE) %>%
        ggplot() +
        geom_point(aes(x=iterations, y=`10-fold cv`), colour=cbp2[2]) +
        geom_line(aes(x=iterations, y=`10-fold cv`), colour=cbp2[2]) +
        geom_ribbon(data=df2, aes(iterations, ymax=qupper, ymin=qlower), fill=cbp2[2], alpha=0.5) +
        geom_vline(xintercept=opt_iter, size=1) +
        xlab("Iteration") +
        ylab("Gaussian loss") +
        ggtitle("Gaussian CV Loss Versus Iterations") +
        theme_bw()

    cat("preparing results...\n")
    res_data <- data.frame(1:niter, cv_nll[1:niter], cv_nll_u[1:niter], cv_nll_l[1:niter])
    names(res_data) <- c("iterations", "cv-mean", "cv_75_quantile", "cv_25_quantile")
    res <- list(
        opt_iter = opt_iter,
        cvplot = p,
        data = res_data
    )
    return(res)

}
