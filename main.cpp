#include <iostream>
#include<math.h>
#include<armadillo>
using namespace std;

const double PI = 3.1415926535987;

float AngWrap (float th);
arma::vec dynRK4_step(arma::vec (*f)(arma::vec, arma::vec), arma::vec& x, arma::vec& u, double dt);

int main()
{
    //const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    
    arma::vec x1 = {1.2, 3.4, 5.6,7.8}; 
    arma::vec x2 = {0.4, 0.3,0.2,0.1};
    //x1+x2;
    cout<<"Hello World!\n";
    //for(int i = 0; i<x1.size();i++){
       cout<<x1+x2;cout<<"\n";
    //}
}

float AngWrap (float th){
    float thwrap;
    thwrap = fmod(th+PI, 2*PI);
    if (thwrap < 0.0) thwrap = thwrap + 2*PI;
    thwrap = thwrap - PI;
    return thwrap;
}

arma::vec dynRK4_step(arma::vec (*f)(arma::vec, arma::vec), arma::vec& x, arma::vec& u, double dt){
    arma::vec k1, k2, k3, k4;
    k1 = f(x,u)*dt;
    k2 = f(x+k1/2, u)*dt;
    k3 = f(x+k2/2, u)*dt;
    k4 = f(x+k3, u)*dt;
    return x + (k1/6)+(k2/3)+(k3/3)+(k4/6);
    
}
//example odeintRK4 (&dynamics, X, U,0.01)