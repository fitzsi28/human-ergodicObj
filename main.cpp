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
    errorcost<CartPend> cost (Q,R,xd,&syst1);
    sac<CartPend,errorcost<CartPend>> sacsys (&syst1,&cost);
    syst1.Ucurr = {-21.0}; 
    syst1.Xcurr = {3.0, 0.0,0.0,0.0};
    myfile<<"time,theta,thetadot,x,xdot,u\n";
    while (syst1.tcurr<5.0){
    myfile<<syst1.tcurr<<",";
    myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<"\n";
    syst1.step();
    }
       
    myfile.close();
    double int1[2]={0.1, 0.2};
    double int2[2] = {0.15,0.18};
    timeInt int3;
    int3 = sacsys.uInterval(int1,int2);
    cout<<int3.start<<", "<<int3.end<<"\n";
     cout<<RK4_step<CartPend>(&syst1, syst1.Xcurr, syst1.Ucurr,0.01)<<"\n";
}

