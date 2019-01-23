#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"cartpend.hpp"
#include"ergodic_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"

double xbound = PI,ybound = 10.;
double xd(double x1, double x2){
  arma::vec lbounds = {{-xbound},{-ybound}};
  arma::vec x = {{x1},{x2}};
  arma::vec Mu = {{0.},{0.}}; Mu=Mu-lbounds;
  arma::mat Sig = {{0.1,0.},{0.,2.}};
return arma::as_scalar(arma::expmat(-0.5*(x-Mu).t()*Sig.i()*(x-Mu))/pow(pow(2*PI,2)*arma::det(Sig),0.5));};
arma::vec unom(double t){
  return arma::zeros(1);};

int main(){
  ofstream myfile;
  myfile.open ("CPergtest.csv");
  CartPend syst1 (1.0,0.1,9.81,2.0,0.01);
  arma::mat R = 0.01*arma::eye(1,1); double q=500.;
  ergodicost<CartPend> cost (q,R,5,0,1,xd,xbound,ybound,1.0,&syst1);
    
  arma::vec umax = {40};
  sac<CartPend,ergodicost<CartPend>> sacsys (&syst1,&cost,0.,2.0,umax,unom);

  arma::vec xwrap;
  syst1.Ucurr = {0.0}; 
  syst1.Xcurr = {3.1, 0.0,0.0,0.0}; 

  myfile<<"time,theta,thetadot,x,xdot,u\n";
  
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
    sacsys.unom_shift(); cost.ckmemory(syst1.Xcurr);
    
    } 
       
    myfile.close();
    ofstream coeff;
    coeff.open("CP_coefficients.csv");
    cost.hk.save(coeff,arma::csv_ascii);
    cost.phik.save(coeff,arma::csv_ascii);
    cost.ckpast.save(coeff,arma::csv_ascii);
    coeff.close();
}

