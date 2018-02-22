#include <iostream>
#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"cartpend.hpp"
#include"error_cost.hpp"

int main()
{   ofstream myfile;
    myfile.open ("test.csv");
    CartPend syst1 (0.1,0.1,9.81,2.0,0.01);
    syst1.Ucurr = {0.0}; 
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
}



