#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"cartpend.hpp"
#include"ergodic_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"

double xd(double x1, double x2){
    arma::vec x = {{x1},{x2}};
    arma::vec Mu = {{PI},{7}};
    arma::mat Sig = {{0.1,0.},{0.,1}};
        return arma::as_scalar(arma::expmat(-0.5*(x-Mu).t()*Sig.i()*(x-Mu))/pow(pow(2*PI,2)*arma::det(Sig),0.5));};
arma::vec unom(double t){
        return arma::zeros(1);};

int main()
{   ofstream myfile;
    myfile.open ("ergtest.csv");
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    arma::mat R = 0.3*arma::eye(1,1); double q=10.;
    double domain[2][2]={{-PI,PI},{-7.,7.}};
    ergodicost<CartPend> cost (q,R,4,xd,domain,&syst1);
 
    arma::vec umax = {20};
    sac<CartPend,ergodicost<CartPend>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
    /*
    arma::mat Q = {
        {200,0.,0.,0.},
        {0., 0.,0.,0.},
        {0.,0.,20.,0.},
        {0.,0.,0.,1.}};
    arma::mat R = 0.3*arma::eye(1,1);
    
 
    arma::vec xwrap;
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {3.1, 0.0,0.0,0.0};
    errorcost<CartPend> cost (Q,R,xd,&syst1);
    sac<CartPend,errorcost<CartPend>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
    arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
       
    myfile<<"time,theta,thetadot,x,xdot,u\n";
 
    while (syst1.tcurr<30.0){
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";//myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";//myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<"\n";
    syst1.step();
    sacsys.SAC_calc();
    syst1.Ucurr = sacsys.ulist.col(0); 
    sacsys.unom_shift();    
    } 
     */  
    myfile.close();
}
