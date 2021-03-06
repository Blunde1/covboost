// Credit to https://gallery.rcpp.org/articles/dmvnorm_arma/
// Modified by Berent Lunde

#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::plugins(openmp)]]


static double const log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
    arma::uword const n = trimat.n_cols;

    for(unsigned j = n; j-- > 0;){
        double tmp(0.);
        for(unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

// [[Rcpp::export(.dmvnorm_arma_mc)]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,
                         arma::rowvec const &mean,
                         arma::mat const &sigma,
                         bool const logd = false,
                         int const cores = 1) {
    using arma::uword;
    omp_set_num_threads(cores);
    uword const n = x.n_rows,
        xdim = x.n_cols;
    arma::vec out(n);

    try{
        arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
        double const rootisum = arma::sum(log(rooti.diag())),
            constants = -(double)xdim/2.0 * log2pi,
            other_terms = rootisum + constants;

        arma::rowvec z;
        #pragma omp parallel for schedule(static) private(z)
        for (uword i = 0; i < n; i++) {
            z = (x.row(i) - mean);
            inplace_tri_mat_mult(z, rooti);
            out(i) = other_terms - 0.5 * arma::dot(z, z);
        }

        if (logd)
            return out;
        return exp(out);
    }  catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rcpp::stop("c++ exception (unknown reason)");
    }

    out.fill(NA_REAL);
    return out;             // not reached


}
