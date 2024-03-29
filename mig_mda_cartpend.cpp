#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"SAC_MDA/cartpend.hpp"
#include"SAC_MDA/error_cost.hpp"
#include"SAC_MDA/SAC.hpp"
#include"SAC_MDA/rk4_int.hpp"
#include"SAC_MDA/MIGMDA.hpp"

arma::vec xd(double t){
        return arma::zeros(4);};
arma::vec unom(double t){
        return arma::zeros(1);};

int main()
{   ofstream myfile;
    myfile.open ("migtest.csv");
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    arma::mat Q = {
        {200,0.,0.,0.},
        {0., 0.,0.,0.},
        {0.,0.,20.,0.},
        {0.,0.,0.,1.}};
    arma::mat R = 0.3*arma::eye(1,1);
    arma::vec umax = {20};
 
    arma::vec xwrap;
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {3.1, 0.0,0.0,0.0};
    errorcost<CartPend> cost (Q,R,xd,&syst1);
    sac<CartPend,errorcost<CartPend>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
    arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
    
    migmda<CartPend,errorcost<CartPend>> demon(&sacsys, false);
    normal_distribution<double> user(0,20);
    default_random_engine generator;
    arma::vec input = {user(generator)};
       
    myfile<<"time,theta,thetadot,x,xdot,u,user\n";
 
    while (syst1.tcurr<30.0){
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";//myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";//myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<","<<input(0)<<"\n";
    syst1.step();
    //sacsys.SAC_calc();
    input = {user(generator)};
    syst1.Ucurr = demon.filter(input);
    //syst1.Ucurr = sacsys.ulist.col(0); 
    sacsys.unom_shift();    
    } 
       
    myfile.close();
}

