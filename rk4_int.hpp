#ifndef RK4_INT_HPP
#define RK4_INT_HPP
#include<armadillo>

arma::vec RK4_step(arma::vec (*f)(arma::vec, arma::vec), arma::vec& x, arma::vec& u, double dt);

#endif