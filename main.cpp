#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"cartpend.hpp"
#include"error_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"

int main()
{   ofstream myfile;
    myfile.open ("test.csv");
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    arma::mat Q = arma::eye(4,4);
    arma::mat R = arma::eye(1,1);
    arma::vec xd = arma::ones(4);
    arma::vec xwrap;
    syst1.Ucurr = {0.0}; 
    syst1.Xcurr = {3.0, 0.0,0.0,0.0};
    errorcost<CartPend> cost (Q,R,xd,&syst1);
    sac<CartPend,errorcost<CartPend>> sacsys (&syst1,&cost);
    //arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
    SACaction u2 = sacsys.SAC_calc(syst1.Xcurr,sacsys.ulist);
    
    myfile<<"time,theta,thetadot,x,xdot,u\n";
 
    while (syst1.tcurr<30.0){
    myfile<<syst1.tcurr<<",";
    //xwrap = syst1.proj_func(syst1.Xcurr);
    myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    //myfile<<xwrap(0)<<","<<xwrap(1)<<",";
    //myfile<<xwrap(2)<<","<<xwrap(3)<<",";
    myfile<<syst1.Ucurr(0)<<"\n";
    syst1.step();
    u2 = sacsys.SAC_calc(syst1.Xcurr,sacsys.ulist);
    //cout<<"\n"<<u2.tau.start<<" "<<syst1.tcurr<<" "<<u2.u;
    syst1.Ucurr = sacsys.ulist(1);
    } 
       
    myfile.close();
    arma::mat xsol; 
    arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
    xsol=sacsys.xforward(unom);
    //cout<<xsol<<"\n";
    //cout<<sacsys.rhoback(xsol,unom);
    //SACaction u2 = sacsys.SAC_calc(syst1.Xcurr,unom);
    //cout<<"\n"<<u2.u;
}

