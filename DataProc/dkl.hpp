#ifndef DKLCOST_HPP
#define DKLCOST_HPP
#include<armadillo>
#include<math.h>
//#include<omp.h>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;

class dklcost {
  double dt;
  double L1=0.5,L2=0.5;
  int K;//number of samples in the search domain
  arma::mat sigma; 
  int MAXINT=240;
  public:
    arma::mat xpast;
    arma::mat domainsamps;
    arma::vec qs_i,ps_i;
    double Q;
    arma::mat R;
    int t_now=0;
    cv::Mat img;
    dklcost(string imageName,int blur, int _K,arma::mat _sigma, double _HIST, double _dt){
      cv::Mat imgtemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
      img = (cv::Scalar::all(255)-imgtemp);
      cv::blur(img,img,cv::Size(blur,blur)); 
      dt=_dt; K = _K;// initialize with dt and  number of samples
      sigma=_sigma;
      MAXINT = _HIST/dt;//cout<<MAXINT<<endl;
      omp_set_dynamic(0); // get rid of dynamic stuff
      omp_set_num_threads(16); // set the number of threads
      xpast.set_size(2,60/dt);//initialize xpast to hold up to a minute of data
      domainsamps.set_size(2,K);
      qs_i.zeros(K);ps_i.set_size(K);
      resample(); //initialize the uniform samples over the domain and the probability of each sample
    };
    double calc_cost (const arma::mat& x);
    void xmemory (const arma::vec&);
    void resample ();
    void qs_disc(const arma::mat& x);
    double phid(const arma::vec& x);
    void resetx ();
};/////////end main class def

double dklcost::calc_cost (const arma::mat& x){
  double J1 = 0.,Jtemp;
  arma::mat xjoined;
  if(t_now<=MAXINT){xjoined = arma::join_rows(xpast.cols(0,t_now),x);}
  else{xjoined = arma::join_rows(xpast.cols(t_now-MAXINT,t_now),x);};
  qs_disc(xjoined);
  J1 = -arma::as_scalar(arma::sum(ps_i%arma::log(qs_i)));
  return J1;}

void dklcost::qs_disc(const arma::mat& x){
  #pragma omp parallel for
  for(int n=0;n<qs_i.n_rows;n++){
    qs_i(n) = 0.;
    for(int j=0;j<x.n_cols;j++){
        arma::vec sj = domainsamps.col(n)-x.col(j);
        qs_i(n)+=dt*exp(-0.5*arma::as_scalar(sj.t()*sigma.i()*sj));
      };
    };
  qs_i=arma::normalise(qs_i,1);//normalise the discrete pdf over the samples
  }
      
void dklcost::resample(){
  //Choose K samples over the domain [[-L1,L1],[-L2,L2]]
  //setup CDF of image for inverse transform sampling
  double epsilon = pow(10,0);//-1);
  double* imgCDF = new double[img.size().width*img.size().height]; 
  double imgNorm = (cv::mean(img)[0]+epsilon)*img.size().width*img.size().height;
    imgCDF[0] = (img.at<uchar>(0,0)+epsilon)/imgNorm; 
    for(int m=0;m<img.size().height;m++){
      for(int n=0;n<img.size().width;n++){
        imgCDF[(m*img.size().height)+n+1] = imgCDF[(m*img.size().height)+n]+(img.at<uchar>(m,n)+epsilon)/imgNorm;
      };
    };//end CDF setup  
  //setup uniform real distribution for inverse transform sampling
  random_device rd; mt19937 eng(rd()); uniform_real_distribution<> distr(0,1);
  
  #pragma omp parallel for
  for(int n=0;n<ps_i.n_rows;n++){
    double k = distr(eng); 
    int i = 0;
    while (imgCDF[i+1]<=k){i++;};
    domainsamps(1,n) = round(i/img.size().width)/img.size().width-0.5;
    domainsamps(0,n) = fmod(i,img.size().width)/img.size().width-0.5;//(i-(i/2200)*2200)-0.5;
    ps_i(n) = phid(domainsamps.col(n));
  };
  ps_i=arma::normalise(ps_i,1);//normalise the discrete pdf over the samples
}

void dklcost::xmemory(const arma::vec& x){
  xpast.col(t_now)= x;
  t_now++;
  //resample();
}
void dklcost::resetx(){
  t_now = 0.;
  //resample();
}

double dklcost::phid(const arma::vec& x){
  double ind1 = (x(1)+L1)*img.size().width; double ind2 = (x(0)+L2)*img.size().width;
  double intensity = img.at<uchar>(round(ind1),round(ind2));
  double totalInt = cv::mean(img)[0]*(L1*2)*(L2*2);
  intensity = intensity/totalInt;
  return intensity;};
#endif