#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;
const int SEARCHRAD = 100;
const int THRESHOLD = 130;

double euclidist (int x, int y, int xt, int yt){
  double d= sqrt(pow(x-xt,2.)+pow(y-yt,2.));
  return d;}


int main()
{ const char* window_name1 = "Original";
  const char* window_name2 = "WithinBounds";
  cv::Mat image;
  string imageName("apple.png");
  cv::Mat imagetemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
  image =(cv::Scalar::all(255)-imagetemp);
  //cv::flip(image,image,-1);
  random_device rd; mt19937 eng(rd());
  uniform_int_distribution<> distr(0,2200-1);
  for(int i=0;i<=100000;i++){
    int x = distr(eng);
    int y = distr(eng);
    int xnear =0,ynear=0;
    double distnear;
    int j = 1; bool found = false;
    while(found==false and j<=SEARCHRAD){
      arma::vec dist = 2200*arma::ones<arma::vec>(8*j);
      arma::mat pixels = arma::zeros<arma::mat>(8*j,2);
      for(int m=-j;m<=j;m++){
        if(y+m>=image.rows or y+m<=0 or x-j>=image.cols or x-j<=0 or x+j>=image.cols or x+j<=0){//do nothing
        }else {
          if(image.at<uchar>(y+m,x-j)>THRESHOLD){found=true;
            dist = arma::shift(dist,1);dist(0) = euclidist(x,y,x-j,y+m);
            pixels = arma::shift(pixels,1,0); pixels(0,0)=x-j; pixels(0,1)=y+m;};
          if(image.at<uchar>(y+m,x+j)>THRESHOLD){found = true;
            dist = arma::shift(dist,1);dist(0) = euclidist(x,y,x+j,y+m);
            pixels = arma::shift(pixels,1,0); pixels(0,0)=x+j; pixels(0,1)=y+m;};
        };
      };
      for(int n=-j+1;n<j;n++){
        if(y-j>=image.rows or y-j<=0 or y+j>=image.rows or y+j<=0 or x+n>=image.cols or x+n<=0){//do nothing
        }else {
          if(image.at<uchar>(y-j,x+n)>THRESHOLD){found=true;
            dist = arma::shift(dist,1);dist(0) = euclidist(x,y,x+n,y-j);
            pixels = arma::shift(pixels,1,0); pixels(0,0)=x+n; pixels(0,1)=y-j;};
          if(image.at<uchar>(y+j,x+n)>THRESHOLD){found = true;
            dist = arma::shift(dist,1);dist(0) = euclidist(x,y,x+n,y+j);
            pixels = arma::shift(pixels,1,0); pixels(0,0)=x+n; pixels(0,1)=y+j;};
        };
      };
      xnear = pixels(dist.index_min(),0);
      ynear = pixels(dist.index_min(),1);
      distnear = euclidist(x,y,xnear,ynear);
      j++;    
    };
    if(found==false){}//imagetemp.at<uchar>(y,x)=0;}
    else{
      if(distnear>50.){imagetemp.at<uchar>(y,x)=0;}
      //cout<<" ("<<x<<","<<y<<") "<<" ("<<xnear<<","<<ynear<<")\n";
      };
  };
 //cv::namedWindow(window_name1,cv::WINDOW_NORMAL);
 cv::namedWindow(window_name2,cv::WINDOW_NORMAL);
 //cv::imshow(window_name1,image);
 cv::imshow(window_name2,imagetemp);
 cv::waitKey(0);
}