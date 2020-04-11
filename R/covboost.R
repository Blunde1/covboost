## Copyright (C) 2020 Berent Lunde
## License: GPL-3

#' Add together two numbers
#'
#' @param x An \code{n x p} numeric matrix or data frame
#' @param learning_rate Scaling the path of elements in the covariance matrix
#' @param niter The number of boosting iterations
#' @return A list of useful boosting information, most importantly the covariance matrix
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @rdname covboost
#' @export
covboost <- function(x, learning_rate=0.01, niter=1000)
{
    # Boosts out a covariance matrix from the Identity matrix
    # for p>>n matrix will become singular. The function notices this and terminates
    # Optimal covariance matrix is likely reached some iterations before this

    # requires: mvtnorm, ggplot2

    # x: data
    # learning_rate: shrink each step
    # niter: max number of iterations

    # see CV version for comments on almost identical code

    cov_heatmap <- function(cov_, title=""){
        # Get upper triangle of the correlation matrix
        get_upper_tri <- function(cormat){
            cormat[lower.tri(cormat)]<- NA
            return(cormat)
        }
        upper_tri <- get_upper_tri(cov_)

        # Melt the correlation matrix
        melted_cormat <- melt(upper_tri, na.rm = TRUE)
        # Heatmap
        ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
            geom_tile(color = "white")+
            scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                 midpoint = 0, limit = c(-1,1), space = "Lab",
                                 name="Pearson\nCorrelation") +
            theme_bw()+
            ggtitle(title) +
            theme(
                legend.justification = c(1, 0),
                legend.position = c(0.6, 0.7),
                legend.direction = "horizontal")#+
        #coord_fixed()
    }


    n <- nrow(x)
    p <- ncol(x)

    nll <- numeric(niter)
    stops <- niter

    B <- diag(p)
    A <- cov(x)

    cat("starting boosting...\n")
    pb <- txtProgressBar(min=0, max=niter, style=3)

    for(i in 1:niter)
    {
        D <- A-B
        ind <- which(abs(D)==max(abs(D)), arr.ind = T)
        B[ind] <- B[ind] + learning_rate*D[ind]
        nll_tmp <- -sum(dmvnorm(x, rep(0,p), B, log = TRUE))

        # checks
        if(!is.finite(nll_tmp)) {
            stops <- i-1
            cat("stopping at iteration: ", i, "\n",
                "updating niter to ", i-1)
            break
        }else{
            #update
            nll[i] <- nll_tmp
        }

        setTxtProgressBar(pb, value=i)
    }
    close(pb)

    p <- cov_heatmap(cov2cor(B), "Boosted Correlation Matrix")

    res <- list(
        cov=B,
        niter=stops,
        plot=p
    )

    return(res)

}
