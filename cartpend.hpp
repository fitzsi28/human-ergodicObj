#ifndef CARTPEND_HPP
#define CARTPEND_HPP
#include<armadillo>
#include"rk4_int.hpp"

const double PI = 3.1415926535987;

class CartPend {
    double m, B, g, h, dt;
    public:
        double tcurr=0.0;
        arma::vec Xcurr, Ucurr;
        CartPend (double, double,double,double,double);
        arma::vec proj_func (const arma::vec& x);
        inline arma::vec f(const arma::vec& x, const arma::vec& u);
        inline arma::mat dfdx(const arma::vec& x, const arma::vec& u);
        arma::vec RK4_int(const arma::vec& x, const arma::vec& u);
        void step(void);
        
};

CartPend::CartPend (double _m, double _B, double _g, double _h, double _dt){
    m = _m; B = _B; g = _g; h=_h;//system parameters
    dt = _dt;//step size
}

arma::vec CartPend::proj_func (const arma::vec& x){
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
    return xwrap;
}
inline arma::vec CartPend::f(const arma::vec& x, const arma::vec& u){
    arma::vec xdot = {x(1),
                      g/h*sin(x(0)) + B*x(1)/(m*h*h)-u(0)*cos(x(0))/h,
                      x(3),
                      u(0)};;
    return xdot;
}; 

inline arma::mat CartPend::dfdx(const arma::vec& x, const arma::vec& u){
    arma::mat A = {
        {0,1,0,0},
        {g/h*cos(x(0))+ u(0)*sin(x(0))/h, B/(m*h*h),0,0},
        {0,0,0,1},
        {0,0,0,0}
    };
    return A;
}; 

arma::vec CartPend::RK4_int(const arma::vec& x, const arma::vec& u){
    arma::vec k1, k2, k3, k4;
    k1 = f(x,u)*dt; 
    k2 = f(x+k1/2, u)*dt; 
    k3 = f(x+k2/2, u)*dt;
    k4 = f(x+k3, u)*dt;
    return (x + (k1/6)+(k2/3)+(k3/3)+(k4/6));
 }

void CartPend::step(){
    Xcurr = RK4_int(Xcurr,Ucurr);
    tcurr = tcurr+dt;
}


#endif