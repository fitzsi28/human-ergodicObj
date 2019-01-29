#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;

#include"doubleint.hpp"
#include"ergodic_cost.hpp"
#include"SAC.hpp"
#include"rk4_int.hpp"

cv::Mat image;
double imgTotal=0.;
double xbound = 456/2,ybound = 599/2;
double phid(double x1, double x2){
  double intensity = image.at<uchar>(round(x2),round(x1));
  //double totalInt = cv::mean(image)[0]*456*599;//cout<<totalInt<<" ";
  intensity = intensity/(0.1*imgTotal);
  return intensity;};

arma::vec unom(double t){
        return arma::zeros(2,1);};

int main()
{   string imageName("Abraham_Lincoln_head_on_shoulders_photo_portraitv3.jpg");
    //string imageName("letter-A.png");
    cv::Mat imagetemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
    image = cv::Scalar::all(255)-imagetemp;
    //double imgTotal = 0.;
    for(int i=0;i<image.size().height;i++){
        for(int j=0;j<image.size().width;j++){
            imgTotal+=(double)image.at<uchar>(i,j);
        }
    }
    xbound = image.size().width/2.; ybound = image.size().height/2.;
    cout<<xbound*2<<" "<<ybound*2<<"\n";
    //cout<<image.size().width<<" "<<image.size().height<<"\n";
    ofstream myfile;
    myfile.open ("DIergtest.csv");
    DoubleInt syst1 (1./60.);
    arma::mat R = 0.1*arma::eye(2,2); double q=1000.;
    arma::vec umax = {200,300};
    double T = 1.0;
    ergodicost<DoubleInt> cost (q,R,10,0,2,phid,xbound,ybound,T,&syst1);
    sac<DoubleInt,ergodicost<DoubleInt>> sacsys (&syst1,&cost,0.,T,umax,unom);
    cout<<q<<"\n";
    arma::vec xwrap;
    syst1.Ucurr = unom(0); 
    syst1.Xcurr = {-100,0.01,-250,0.01};
    
    //arma::mat unom = arma::zeros<arma::mat>(1,sacsys.T_index);
       
    myfile<<"time,x,xdot,y,ydot,ux,uy\n";
 
    while (syst1.tcurr<60.0){
    if(fmod(syst1.tcurr,5)<syst1.dt)cout<<"Time: "<<syst1.tcurr<<"\n";
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";//myfile<<syst1.Xcurr(0)<<","<<syst1.Xcurr(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";//myfile<<syst1.Xcurr(2)<<","<<syst1.Xcurr(3)<<",";
    myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<"\n";
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

