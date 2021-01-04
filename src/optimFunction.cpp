#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <math.h>       /* sqrt */
#include "roptim.h"
// [[Rcpp::depends(roptim)]]


using namespace Rcpp;
using namespace std;
using namespace roptim;
//' bfunction
//' @param X high dimensional matrix
//' @param theta unknown regression vector
//' @param family three family cases: gaussian, poisson and binomial
//' @export
// [[Rcpp::export]]
arma::vec bFunction (arma::mat X, arma::vec theta,const std::string family){
  arma::vec result;
  if (family == "gaussian"){
    result = 0.5 * (X * theta) % (X * theta);
  } else if (family == "poisson"){
    result = exp(X * theta);
  } else {
    result = log(1 + exp(X * theta));
  }
  return result;
}
//' mean function
//' @param X high dimensional matrix
//' @param theta unknown regression vector
//' @param family three family cases: gaussian, poisson and binomial
//' @export
// [[Rcpp::export]]
arma::vec MeanFunction (arma::mat X, arma::vec theta,const std::string family){
  arma::vec result;
  if (family == "gaussian"){
    result = X * theta;
  } else if (family == "poisson"){
    result = exp(X * theta);
  } else {
    result = exp(X * theta) / ( 1 + exp(X * theta));
  }
  return result;
}
//' mean prime function
//' @param X high dimensional matrix
//' @param theta unknown regression vector
//' @param family three family cases: gaussian, poisson and binomial
//' @export
// [[Rcpp::export]]
arma::vec MeanPrime(arma::mat X, arma::vec theta,const std::string family){
  arma::vec result;
  if (family == "gaussian"){
    arma::vec A(X.n_cols, arma::fill::ones);
    result = A;
  } else if (family == "poisson"){
    result = exp(X * theta);
  } else {
    result = exp(X * theta) / pow((1 +
      exp(X * theta)) , 2);
  }
  return result;
}
//' objective function
//' @param theta unknown regression vector
//' @param family three family cases: gaussian, poisson and binomial
//' @param y where y = Xβ + ε
//' @param X high dimensional matrix
//' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
//' @export
// [[Rcpp::export]]
double ObjectiveFunction(arma::vec theta,const std::string family,arma::vec y,arma::mat X,double trex_c) {
  if(family == "gaussian")
  { arma::vec loss = y - X*theta;
    arma::vec derivative = X.t() * (y - X *theta);
    double result = (pow(norm(loss, 2) , 2)- pow(norm(y, 2) , 2)) / (trex_c * norm(derivative, 2)) + norm(theta, 2);
    return result;
  }
  else
  { double loss = sum(y % (X * theta) - bFunction(X, theta,family));
    arma::vec tune_vector = ((y - MeanFunction(X, theta,family)).t() * X).t();
    double regularization = norm(theta, 2);
    double result = -loss /  (trex_c * norm(tune_vector,2 )) + regularization;
    
    return result;
  }
}
//' Gradient of objective function
//' @param theta unknown regression vector
//' @param family three family cases: gaussian, poisson and binomial
//' @param y where y = Xβ + ε
//' @param X high dimensional matrix
//' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
//' @export
// [[Rcpp::export]]
arma::vec Gradient (arma::vec theta,const std::string family,arma::vec y,arma::mat X,double trex_c) {
  if(family == "gaussian") {
    arma::vec derivative = X.t() * (y - X * theta);
    arma::vec loss = (y - X * theta);
    arma::vec result = (pow(norm(loss, 2) , 2)- pow(norm(y, 2) , 2)) * X.t() * X * derivative /
      (trex_c * pow(norm(derivative, 2) , 3)) -
        2 * derivative / (trex_c * norm(derivative,2 )) +
        theta / norm(theta, 2);
    return result;
  } else {
    arma::vec meanprime = MeanPrime(X, theta,family);
    arma::mat MP(X.n_rows, X.n_cols, arma::fill::zeros);
    for (int i = 0; i < X.n_cols; ++i) {
      MP.col(i)=meanprime;}
    double loss = sum(y % (X * theta) - bFunction(X, theta,family));
    arma::vec tune_vector = ((y - MeanFunction(X, theta,family)).t() * X).t();
    arma::vec result = -tune_vector / (trex_c * norm(tune_vector,2 )) -
      (0.5 * trex_c * 2 * loss * (tune_vector).t() * X.t() * (X % MP) / (pow(trex_c , 2) * pow(norm(tune_vector, 2) , 3))).t() +
      theta / norm(theta, 2);
    return result;
  }
}
//' GradientLs function
//' @param theta unknown regression vector
//' @param X high dimensional matrix
//' @param y where y = Xβ + ε
//' @param family three family cases: gaussian, poisson and binomial
//' @export
// [[Rcpp::export]]
arma::vec GradientLs (arma::vec theta,arma::mat X,arma::vec y,const std::string family) {
  if(family=="gaussian") {
    arma::vec derivative = X.t() * (y - X * theta);
    arma::vec result = -2 * derivative;
    return result;
  } else {
    arma::vec tune_vector = ((y - MeanFunction(X, theta,family)).t() * X).t();
    arma::vec result = -tune_vector;
    return result;
  }
}
class ObLs : public Functor {
public:
  ObLs(const arma::mat &X,const arma::vec &y,const std::string & family) : X_(X),y_(y),family_(family) {}
  
  double operator()(const arma::vec &theta) override {
    if(family_=="gaussian") {
      arma::vec loss =y_ - X_ * theta;
      double result =  pow(norm(loss, 2) , 2)- pow(norm(y_, 2) , 2);
      return result;
    } else {
      double loss = sum(y_ % (X_ * theta) - bFunction(X_, theta,family_));
      double result = -loss;
      return result;
    }
  }
  void Gradient(const arma::vec &theta, arma::vec &grad) override {
    grad = theta;
    if(family_=="gaussian") {
      arma::vec derivative = X_.t() * (y_ - X_ * theta);
      arma::vec result = -2 * derivative;
      grad = result;
    } else {
      arma::vec tune_vector = ((y_ - MeanFunction(X_, theta,family_)).t() * X_).t();
      arma::vec result = -tune_vector;
      grad= result;
    }
  }
private:
  arma::mat X_;
  arma::vec y_;
  std::string family_;
  
};
//' find optimized solution for ObjectLs function
//' @param theta vector store our optimized result
//' @param X high dimensional matrix
//' @param y where y = Xβ + ε
//' @param family three family cases: gaussian, poisson and binomial
//' @export
// [[Rcpp::export]]
arma::vec optim_ObLs(arma::vec theta,arma::mat X,arma::vec y,const std::string family){
  ObLs dist(X,y,family);
  Roptim<ObLs> opt("CG");
  opt.control.type = 1;
  opt.minimize(dist, theta);
  return theta;
  
}


class ObRidge : public Functor {
public:
  ObRidge(const arma::mat &X,const arma::vec &y,const std::string & family,const double &r) : X_(X),y_(y),family_(family),r_(r) {}
  
  double operator()(const arma::vec &theta) override {
    if(family_=="gaussian") {
      arma::vec loss =y_ - X_ * theta;
      double result =  pow(norm(loss, 2) , 2)- pow(norm(y_, 2) , 2) +r_*pow(norm(theta, 2) , 2);
      return result;
    } else {
      double loss = sum(y_ % (X_ * theta) - bFunction(X_, theta,family_));
      double result = -loss+r_*pow(norm(theta, 2) , 2);
      return result;
    }
  }
  void Gradient(const arma::vec &theta, arma::vec &grad) override {
    grad = theta;
    if(family_=="gaussian") {
      arma::vec derivative = X_.t() * (y_ - X_ * theta);
      arma::vec result = -2 * derivative + r_ * 2 *norm(theta, 2);
      grad = result;
    } else {
      arma::vec tune_vector = ((y_ - MeanFunction(X_, theta,family_)).t() * X_).t();
      arma::vec result = -tune_vector + r_ * 2 *norm(theta, 2);
      grad= result;
    }
  }
  
private:
  arma::mat X_;
  arma::vec y_;
  std::string family_;
  double r_;
};

//' find optimized solution for ObjectRidge function
//' @param theta vector store our optimized result
//' @param X high dimensional matrix
//' @param y where y = Xβ + ε
//' @param family three family cases: gaussian, poisson and binomial
//' @param r tuning.parameter
//' @export
// [[Rcpp::export]]
arma::vec optim_Ridge(arma::vec theta,arma::mat X,arma::vec y,const std::string family,double r){
  ObRidge dist3(X,y,family,r);
  Roptim<ObRidge> opt3("CG");
  opt3.control.type = 1;
  opt3.minimize(dist3, theta);
  return theta;
  
}



class ObFn : public Functor {
public:
  ObFn(const arma::mat &X,const arma::vec &y,const std::string & family,const double &trex_c) : X_(X),y_(y),family_(family),
  trex_c_(trex_c){}
  
  double operator()(const arma::vec &theta) override {
    if(family_ == "gaussian")
    { arma::vec loss = y_ - X_*theta;
      arma::vec derivative = X_.t() * (y_ - X_ *theta);
      double result = (pow(norm(loss, 2) , 2)- pow(norm(y_, 2) , 2)) / (trex_c_ * norm(derivative, 2)) + norm(theta, 2);
      return result;
    }
    else
    { double loss = sum(y_ % (X_ * theta) - bFunction(X_, theta,family_));
      arma::vec tune_vector = ((y_ - MeanFunction(X_, theta,family_)).t() * X_).t();
      double regularization = norm(theta, 2);
      double result = -loss /  (trex_c_ * norm(tune_vector,2 )) + regularization;
      return result;
    }
  }
  void Gradient(const arma::vec &theta, arma::vec &grad) override {
    grad = theta;
    if(family_ == "gaussian") {
      arma::vec derivative = X_.t() * (y_ - X_ * theta);
      arma::vec loss = (y_ - X_ * theta);
      arma::vec result = (pow(norm(loss, 2) , 2)- pow(norm(y_, 2) , 2)) * X_.t() * X_ * derivative /
        (trex_c_ * pow(norm(derivative, 2) , 3)) -
          2 * derivative / (trex_c_ * norm(derivative,2 )) +
          theta / norm(theta, 2);
      grad= result;
    } else {
      arma::vec meanprime = MeanPrime(X_, theta,family_);
      arma::mat MP(X_.n_rows, X_.n_cols, arma::fill::zeros);
      for (int i = 0; i < X_.n_cols; ++i) {
        MP.col(i)=meanprime;}
      double loss = sum(y_ % (X_ * theta) - bFunction(X_, theta,family_));
      arma::vec tune_vector = ((y_ - MeanFunction(X_, theta,family_)).t() * X_).t();
      arma::vec result = -tune_vector / (trex_c_ * norm(tune_vector,2 )) -
        (0.5 * trex_c_ * 2 * loss * (tune_vector).t() * X_.t() * (X_ % MP) / (pow(trex_c_ , 2) * pow(norm(tune_vector, 2) , 3))).t() +
        theta / norm(theta, 2);
      grad= result;
    }
  }
private:
  arma::mat X_;
  arma::vec y_;
  std::string family_;
  double trex_c_;
  
};
//' find optimized solution for Objective function
//' @param theta vector store our optimized result
//' @param X high dimensional matrix
//' @param y where y = Xβ + ε
//' @param family three family cases: gaussian, poisson and binomial
//' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
//' @export
// [[Rcpp::export]]
arma::vec optim_ObFn(arma::vec theta,arma::mat X,arma::vec y,const std::string family,double trex_c){
  ObFn dist2(X,y,family,trex_c);
  Roptim<ObFn> opt2("CG");
  opt2.control.type = 1;
  opt2.minimize(dist2, theta);
  
  return theta;
  
}
