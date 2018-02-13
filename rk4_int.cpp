#include"rk4_int.hpp"

arma::vec RK4_step(arma::vec (*f)(arma::vec, arma::vec), arma::vec& x, arma::vec& u, double dt){
    arma::vec k1, k2, k3, k4;
    k1 = f(x,u)*dt; 
    k2 = f(x+k1/2, u)*dt; 
    k3 = f(x+k2/2, u)*dt;
    k4 = f(x+k3, u)*dt;
    return x + (k1/6)+(k2/3)+(k3/3)+(k4/6);
    
};//example odeintRK4 (&dynamics, X, U,0.01)
