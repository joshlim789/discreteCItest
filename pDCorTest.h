#ifndef PDCORTEST_H
#define PDCORTEST_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

StringVector matrix_to_string (StringMatrix sep_vectors);

double get_G2_one(StringVector A, StringVector B, int tot_Au_size, int tot_Bu_size);

double get_G2_all(StringVector A, StringVector B, StringVector S);

List condInttestdis(StringMatrix df, const size_t &i,const size_t &j,
                     const arma::uvec &k, const double &signif_level);

#endif