## Copyright (C) 2020 Berent Lunde
## License: GPL-3

#' Boosted Covariance Matrix Estimation
#'
#' @param x An \code{n x p} numeric matrix or data frame
#' @param shrinkage Shrinkage defining the trade-off between the Identity and variance matrices
#' @param learning_rate Scaling the path of elements in the covariance matrix
#' @param niter The number of boosting iterations
#' @return A list of useful boosting information, most importantly the covariance matrix
#' @examples
#' \dontrun{
#' ## Generate data with random correlation and estimate
#' p <- 50
#' n <- 50
#' x <- matrix(nrow=n, ncol=p)
#' x[,1] <- rnorm(n)
#' for(j in 2:p){
#'     rho <- runif(1,-1,1)*rbinom(1,1,0.5)
#'     col_dependence <- sample(j-1,1)
#'     x[,j] <- rnorm(n, mean=rho*x[,col_dependence], sd=(1-rho^2))
#' }
#' lrn_rate <- 0.2
#' cov_cv <- covboost_cv(x, learning_rate=lrn_rate)
#' cov_cv$cvplot
#' cov_cv$opt_iter
#' cov_cv$opt_shrinkage
#' sigma <- covboost(x, shrinkage=cov_cv$opt_shrinkage, niter=cov_cv$opt_iter, learning_rate=lrn_rate)
#' sigma$cov
#' sigma$plot
#' }
#' @importFrom stats cov cov2cor quantile
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @rdname covboost
#' @export
covboost <- function(x, shrinkage=0.1, learning_rate=0.01, niter=1000)
{
    # Boosts out a covariance matrix from the Identity matrix
    # for p>>n matrix will become singular. The function notices this and terminates
    # Optimal covariance matrix is likely reached some iterations before this

    # requires: mvtnorm, ggplot2

    # x: data
    # learning_rate: shrink each step
    # niter: max number of iterations

    # see CV version for comments on almost identical code

    `10-fold cv` <- cv <- Var1 <- Var2 <- iterations <- qlower <- qupper <- value <- NULL

    n <- nrow(x)
    p <- ncol(x)

    nll <- numeric(niter)
    stops <- niter+1

    # Matrices
    I <- B_rho <- diag(p)
    dA <- cov(x) * shrinkage
    dV <- diag(diag(dA))

    # Standardization
    dA_rho <- cov2cor(dA)

    cat("starting boosting...\n")
    pb <- txtProgressBar(min=0, max=niter+1, style=3)

    if(niter>1){
        for(i in 2:(niter+1))
        {
            # max diff
            dD_rho <- dA_rho - B_rho
            ind <- which(abs(dD_rho)==max(abs(dD_rho)), arr.ind=T)

            # update
            B_rho[ind] <- B_rho[ind] + learning_rate*dD_rho[ind]

            # Does not need check -- I think...!

            setTxtProgressBar(pb, value=i)
        }
        close(pb)
    }

    # transform to final matrix
    B <- sqrt(dV) %*% B_rho %*% sqrt(dV)

    if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
        .cov_heatmap <- function(cov_, title=""){
            # Get upper triangle of the correlation matrix
            get_upper_tri <- function(cormat){
                cormat[lower.tri(cormat)]<- NA
                return(cormat)
            }
            upper_tri <- get_upper_tri(cov_)

            # Melt the correlation matrix
            melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
            # Heatmap
            ggplot2::ggplot(data = melted_cormat, ggplot2::aes(Var2, Var1, fill = value))+
                ggplot2::geom_tile(color = "white")+
                ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                     midpoint = 0, limit = c(-1,1), space = "Lab",
                                     name="Pearson\nCorrelation") +
                ggplot2::theme_bw()+
                ggplot2::ggtitle(title) +
                ggplot2::theme(
                    legend.justification = c(1, 0),
                    legend.position = c(0.6, 0.7),
                    legend.direction = "horizontal")#+
            #coord_fixed()
        }

        p <- .cov_heatmap(cov2cor(B), "Boosted Correlation Matrix")

    } else {
        p <- cat("Install ggplo2 and reshape2 packages to get cov boosting \n")
    }


    res <- list(
        cov=B,
        cor=B_rho,
        niter=stops,
        plot=p
    )

    return(res)

}
