#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;

#include"SAC_MDA/doubleint.hpp"
#include"SAC_MDA/ergodic_cost.hpp"
#include"SAC_MDA/SAC.hpp"
#include"SAC_MDA/rk4_int.hpp"

cv::Mat image;
double imgTotal=0.;
double xbound = 0.5,ybound = 0.5;//2200/2;
double phid(double x1, double x2){
  double ind1 = x2*2200.; double ind2 = x1*2200.;
  double intensity = image.at<uchar>(round(ind1),round(ind2));
  double totalInt = cv::mean(image)[0]*(xbound*2)*(ybound*2);//cout<<totalInt<<" ";
  intensity = intensity/totalInt;//(255*7);
  return intensity;};

arma::vec unom(double t){
        return arma::zeros(2,1);};

int main()
{   string imageName("lincoln2.png");
    //string imageName("gauss.png");
    //string imageName("apple.png");
    cv::Mat imagetemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
    image = (cv::Scalar::all(255)-imagetemp);cout<<cv::mean(image)[0]<<"\n";
    cv::flip(image,image,-1);
    cout<<image.size().width<<" "<<image.size().height<<"\n"; 
    //double imgTotal = 0.;
    //xbound = image.size().width/2.; ybound = image.size().height/2.;
    cout<<xbound*2<<" "<<ybound*2<<"\n";
    ofstream myfile;
    myfile.open ("DIergtest.csv");
    DoubleInt syst1 (1./60.);
    arma::mat R = 0.01*arma::eye(2,2); double q=2000.;
    arma::vec umax = {40,40};
    double T = 1.0;
    ergodicost<DoubleInt> cost (q,R,10,0,2,phid,xbound,ybound,T,&syst1);
    sac<DoubleInt,ergodicost<DoubleInt>> sacsys (&syst1,&cost,0.,T,umax,unom);
    arma::vec xwrap;
    syst1.Ucurr = unom(0); 
    random_device rd; mt19937 eng(rd());
    uniform_real_distribution<> distr(-0.4,0.4);
    syst1.Xcurr = {distr(eng),distr(eng),distr(eng),distr(eng)};
    cout<<syst1.Xcurr<<"\n";
    //arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
       
    myfile<<"time,x,xdot,y,ydot,ux,uy,ergcost\n";
 
    while (syst1.tcurr<30.){
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";//myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";//myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<",";
    myfile<<cost.calc_cost(syst1.Xcurr,syst1.Ucurr)<<"\n";
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

