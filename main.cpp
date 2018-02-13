#include <iostream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"rk4_int.hpp"
const double PI = 3.1415926535987;

float AngWrap (float th);
arma::vec sysdyn(arma::vec x, arma::vec u);

int main()
{
    //const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    
    arma::vec u = {0.0}; 
    arma::vec x1 = {3.1, 0.0,0.0,0.0};
    arma::vec x2;
    x2=RK4_step(&sysdyn, x1, u,0.01);
    cout<<"Hello World!\n";
    cout<<sysdyn(x1,u);cout<<"\n";   
    cout<<x2;cout<<"\n";
    
}

float AngWrap (float th){
    float thwrap;
    thwrap = fmod(th+PI, 2*PI);
    if (thwrap < 0.0) thwrap = thwrap + 2*PI;
    thwrap = thwrap - PI;
    return thwrap;
}

arma::vec sysdyn(arma::vec x, arma::vec u){
    const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    arma::vec xdot = {x(1),
                      g/(h*sin(x(0))) + B*x(1)/(m*h*h*sin(x(0))*sin(x(0)))-x(1)*x(1)/tan(x(0)),
                      x(3),
                      g/tan(x(0))+B*x(1)/(m*h*tan(x(0))*sin(x(0)))-h*x(1)*x(1)/sin(x(0))};
    return xdot;
} 

 /*arma::vec sysdyn(arma::vec& x, arma::vec& u){
    const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    arma::vec xdot = {x(1),
                      g/h*sin(x(0)) + B*x(1)/(m*h*h)-u(0)*cos(x(0))/h,
                      x(3),
                      u(0)};
     return xdot;
} */