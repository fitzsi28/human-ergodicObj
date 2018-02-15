#include <iostream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"rk4_int.hpp"
#include"cartpend.hpp"
const double PI = 3.1415926535987;

arma::vec proj_func (arma::vec& x);
arma::vec sysdyn(const arma::vec& x, const arma::vec& u);

int main()
{
    //const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    CartPend syst1 (0.1,0.01,9.81,2.0);
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {3.6, 0.0,0.0,0.0};
    arma::vec x2;
    x2=RK4_step<CartPend>(syst1, syst1.Xcurr, syst1.Ucurr,0.01);
    //syst1.rk_step()
    cout<<"Hello World!\n";
    cout<<x2;cout<<"\n";   
    cout<<proj_func(syst1.Xcurr);cout<<"\n";
    
}

arma::vec proj_func (arma::vec& x){
    arma::vec xwrap=x;
    xwrap(0) = fmod(x(0)+PI, 2*PI);
    if (xwrap(0) < 0.0) xwrap(0) = xwrap(0) + 2*PI;
    xwrap(0) = xwrap(0) - PI;
    return xwrap;
}

