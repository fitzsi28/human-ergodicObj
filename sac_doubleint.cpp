#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"doubleint.hpp"
#include"error_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"
arma::vec xd(double t){
        return {{-2},{0},{-2},{0}};};
arma::vec unom(double t){
        return arma::zeros(2,1);};

int main()
{   ofstream myfile;
    myfile.open ("DItest.csv");
    DoubleInt syst1 (0.01);
    arma::mat Q = {
        {100,0.,0.,0.},
        {0., 10.,0.,0.},
        {0.,0.,100.,0.},
        {0.,0.,0.,10.}};
    arma::mat R = 0.3*arma::eye(2,2);
    arma::vec umax = {10,10};
 
    arma::vec xwrap;
    syst1.Ucurr = {{0.0},{0.0}}; 
    syst1.Xcurr = {1.1,-0.1,1,0.1};
    errorcost<DoubleInt> cost (Q,R,xd,&syst1);
    sac<DoubleInt,errorcost<DoubleInt>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
    //arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
       
    myfile<<"time,x,xdot,y,ydot,u\n";
 
    while (syst1.tcurr<30.0){
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
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
      
    myfile.close();
    
}

