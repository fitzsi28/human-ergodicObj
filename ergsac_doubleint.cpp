#include <iostream>

#include <fstream>
#include<math.h>
#include<armadillo>
using namespace std;

#include"doubleint.hpp"
#include"ergodic_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"

double xbound = 5,ybound = 5;
double phid(double x1, double x2){
  arma::vec lbounds = {{-xbound},{-ybound}};
  arma::vec x = {{x1},{x2}};
  arma::vec Mu= {{-2.},{1.}}; Mu=Mu-lbounds;
  arma::vec Mu2 = {{3.},{-1.}}; Mu2 = Mu2-lbounds;
  arma::mat Sig = {{0.5,0.},{0.,0.5}};
  double dist1 = arma::as_scalar(arma::expmat(-0.5*(x-Mu).t()*Sig.i()*(x-Mu))/pow(pow(2*PI,2)*arma::det(Sig),0.5));
  double dist2 = arma::as_scalar(arma::expmat(-0.5*(x-Mu2).t()*Sig.i()*(x-Mu2))/pow(pow(2*PI,2)*arma::det(Sig),0.5));
return dist1+dist2;};
arma::vec unom(double t){
        return arma::zeros(2,1);};

int main()
{   ofstream myfile;
    myfile.open ("DIergtest.csv");
    DoubleInt syst1 (1./60.);
    arma::mat Q = {
        {1.,0.,0.,0.},
        {0.,0.,0.,0.},
        {0.,0.,1.,0.},
        {0.,0.,0.,0.}};
    arma::mat R = 0.3*arma::eye(2,2); double q=1000.;
    arma::vec umax = {5,5};
    ergodicost<DoubleInt> cost (q,Q,R,10,0,2,phid,xbound,ybound,1.0,&syst1);
    sac<DoubleInt,ergodicost<DoubleInt>> sacsys (&syst1,&cost,0.,1.0,umax,unom);
 
    arma::vec xwrap;
    syst1.Ucurr = unom(0); 
    syst1.Xcurr = {0,-0.1,0,0.1};
    
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
    sacsys.unom_shift(); cost.ckmemory(syst1.Xcurr);  
    } 
      
    myfile.close();
 ofstream coeff;
 coeff.open("DI_coefficients.csv");
 cost.hk.save(coeff,arma::csv_ascii);
 cost.phik.save(coeff,arma::csv_ascii);
 cost.ckpast.save(coeff,arma::csv_ascii);
 coeff.close();
    
}

