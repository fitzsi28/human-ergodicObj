#ifndef SINGLEINT_HPP
#define SINGLEINT_HPP
#include<armadillo>
#include"rk4_int.hpp"

const double PI = 3.1415926535987;

class SingleInt {
    public:
        double dt,B=0.01;
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        arma::mat xdlist;
        SingleInt (double);
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        inline arma::mat hx(const arma::vec& x);
        void step(void);
        
        
};

SingleInt::SingleInt (double _dt){
   dt = _dt;//step size
}

arma::vec SingleInt::proj_func (const arma::vec& x){
    return x;
}
inline arma::vec SingleInt::f(const arma::vec& x, const arma::vec& u){
    arma::vec xdot = {u(0),
                      u(1)}; 
    return xdot;
}; 

inline arma::mat SingleInt::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {
        {0,0},
        {0,0}
    };
    return A;
}; 

inline arma::mat SingleInt::hx(const arma::vec& x){
    arma::mat H = {
        {1,0},
        {0,1}
    };
    return H;
}; 

void SingleInt::step(){
    Xcurr = RK4_step(this,Xcurr,Ucurr,dt);
    tcurr = tcurr+dt;
};


#endif