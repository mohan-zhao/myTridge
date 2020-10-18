#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <math.h>       /* sqrt */



using namespace Rcpp;
using namespace std;
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat mvrnormArma(const int n, arma::vec mu, const int p, const double rho) {
  
  arma::mat sigma(p, p, arma::fill::zeros);
  
  for (int i = 0; i < sigma.n_rows; ++i) {
    for (int j = 0; j < sigma.n_cols; ++j) {
      sigma(i,j) = pow(rho, abs((i + 1) - (j + 1)));
    }
  }
  
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  /*return arma::normalise(X,2,0);*/
  
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::vec rpois_rcpp( arma::vec &lambda) {
  int n= lambda.n_elem;
  unsigned int lambda_i = 0;
  IntegerVector sim(n);
  for (unsigned int i = 0; i < n; i++) {
    sim[i] = R::rpois(lambda[lambda_i]);
    // update lambda_i to match next realized value with correct mean
    lambda_i++;
  }
  return  as<arma::vec>(sim);
}

// [[Rcpp::export]]
arma::vec cpprbinom(int n, double size, arma::vec prob) {
  NumericVector v = no_init(n);
  std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(size, p); });
  return as<arma::vec>(v);
}

// [[Rcpp::export]]
List genDataList(const int n, const arma::vec& mu, int p, double rho,
                 arma::vec& beta, const double SNR,const std::string Test_case) {
  
  arma::mat U, V, data, normData, Projection;
  arma::vec s, y, means, noise;
  data = mvrnormArma(n, mu, p, rho);
  normData = arma::normalise(data,2,0);
  arma::svd_econ(U,s,V,normData,"right");
  Projection = V * trans(V);
  beta = Projection * beta;
  
  if(Test_case == "gaussian")
  {
    means=normData * beta;
    y = means + arma::randn(n) * sqrt(arma::var(means) / SNR);
  }
  else if (Test_case == "poisson")
  {
    means=exp(normData * beta);
    y = rpois_rcpp(means);
  }
  else
  {
    means=exp(normData * beta)/(1 + exp(normData * beta));
    y = cpprbinom(n,1,means);
  }
  
  
  List ret;
  ret["normData"] = normData;
  ret["beta"] = beta;
  ret["y"] = y;
  ret["Test.case"]=Test_case;
  return ret;
  
}
