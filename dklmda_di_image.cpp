#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;

#include"SAC_MDA/doubleint.hpp"
#include"SAC_MDA/dklimg_cost.hpp"
#include"SAC_MDA/SAC.hpp"
#include"SAC_MDA/rk4_int.hpp"
#include"SAC_MDA/MIGMDA.hpp"


arma::vec unom(double t){
        return arma::zeros(2,1);};

int main()
{   double xbound = 0.5,ybound = 0.5;
    cv::Mat image;
    //string imageName("lincoln2.png");
    //string imageName("gauss.png");
    //string imageName("apple.png");
    string imageName("house.png");
    cv::Mat imagetemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
    image = (cv::Scalar::all(255)-imagetemp);
    cv::flip(image,image,0);
    cv::blur(image,image,cv::Size(50,50));//apple&house 50, carv4&umbrella 100,banana 0
    
    ofstream myfile;
    myfile.open ("DIdkltest.csv");
    DoubleInt syst1 (1./60.);
    syst1.Ucurr = unom(0); 
    random_device rd; mt19937 eng(rd());
    uniform_real_distribution<> distr(-0.4,0.4);
    //syst1.Xcurr = {-0.2,0.0,0.1,0.0};
    syst1.Xcurr = {distr(eng),distr(eng),distr(eng),distr(eng)};//must be initialized before instantiating cost
    arma::mat R = 0.0001*arma::eye(2,2); double q=1.;
    arma::vec umax = {40.0,40.0};
    double T = 0.7;//0.6;
    arma::mat SIGMA = 0.01*arma::eye(2,2);
    dklcost<DoubleInt> cost (q,R,65,SIGMA,0,2,image,xbound,ybound,T,4.0,&syst1);
    sac<DoubleInt,dklcost<DoubleInt>> sacsys (&syst1,&cost,0.,T,umax,unom);
    migmda<DoubleInt,dklcost<DoubleInt>> demon(&sacsys, false);
    arma::vec xwrap;
    uniform_real_distribution<> user(-10,10);
    default_random_engine generator;
    arma::vec input = {user(generator),user(generator)};
    
    myfile<<"time,x,xdot,y,ydot,ux,uy,dklcost\n";
    double start_time = omp_get_wtime();
    while (syst1.tcurr<10.0){
    
    cost.xmemory(syst1.Xcurr);
    //cout <<"resamp time: "<< 1000 * (omp_get_wtime() - start_time)<<endl;
    if(fmod(syst1.tcurr,2)<syst1.dt){
        cout<<"Time: "<<syst1.tcurr<<"\n";
        cout<<"Time Elapsed: "<<omp_get_wtime()-start_time<<endl;
        start_time = omp_get_wtime();
    }
    myfile<<syst1.tcurr<<",";
    xwrap = syst1.proj_func(syst1.Xcurr); 
    myfile<<xwrap(0)<<","<<xwrap(1)<<",";
    myfile<<xwrap(2)<<","<<xwrap(3)<<",";
    myfile<<syst1.Ucurr(0)<<","<<syst1.Ucurr(1)<<",";
    //myfile<<cost.calc_cost(syst1.Xcurr,syst1.Ucurr);
    myfile<<"\n";
    syst1.step();
    //sacsys.SAC_calc();
    input = {user(generator),user(generator)};
    syst1.Ucurr = demon.filter(input); 
    sacsys.unom_shift();  
     
    } 
      
 myfile.close();
 ofstream samples;
 samples.open("Domain_samples.csv");
 cost.domainsamps.save(samples,arma::csv_ascii);
 arma::mat temp = cost.ps_i.t();
 temp.save(samples,arma::csv_ascii);
 samples.close();
 
 //cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );// Create a window for display.
 //cv::imshow( "Display window", image );                   // Show our image inside it.
 //cv::waitKey(0); 
}

