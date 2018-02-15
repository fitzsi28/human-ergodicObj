#ifndef CARTPEND_HPP
#define CARTPEND_HPP
#include<armadillo>


class CartPend {
    float m, B, g, h;
    public:
        arma::vec Xcurr, Ucurr;
        CartPend (float, float,float,float);
        arma::vec f(const arma::vec& x, const arma::vec& u);
        
};

CartPend::CartPend (float a, float b, float c, float d){
    m = a; B = b; g = c; h=d;
}

arma::vec CartPend::f(const arma::vec& x, const arma::vec& u){
    //const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    arma::vec xdot = {x(1),
                      g/(h*sin(x(0))) + B*x(1)/(m*h*h*sin(x(0))*sin(x(0)))-x(1)*x(1)/tan(x(0)),
                      x(3),
                      g/tan(x(0))+B*x(1)/(m*h*tan(x(0))*sin(x(0)))-h*x(1)*x(1)/sin(x(0))};
    return xdot;
} 


#endif