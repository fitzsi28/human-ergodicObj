#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
//#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;
const int SEARCHRAD = 50;
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
  image = (cv::Scalar::all(255)-imagetemp);//cout<<cv::mean(image)[0]<<"\n";
  //cv::flip(image,image,-1);
  random_device rd; mt19937 eng(rd());
  uniform_int_distribution<> distr(0,2200);//-SEARCHRAD);
  //cout<<"\n"<<(double)imagetemp.at<uchar>(2200,distr(eng))<<"\n";
  //cout<<(double)imagetemp.at<uchar>(0,10)<<"\n";
  for(int i=0;i<=100000;i++){
    int x = distr(eng);
    int y = distr(eng);
    int j = 1;
    bool found = false;
    while(found==false and j<=SEARCHRAD){//for(int j=1;j<=SEARCHRAD;j++){
      for(int m=-j;m<=j;m++){
        if(y+m>=image.rows or y+m<=0 or x-j>=image.cols or x-j<=0 or x+j>=image.cols or x+j<=0){//do nothing
        }else {
          if(image.at<uchar>(y+m,x-j)>THRESHOLD){imagetemp.at<uchar>(y,x)=0;found=true;};
          if(image.at<uchar>(y+m,x+j)>THRESHOLD){imagetemp.at<uchar>(y,x)=0;found = true;};
        };
      };
      for(int n=-j+1;n<j;n++){
        if(y-j>=image.rows or y-j<=0 or y+j>=image.rows or y+j<=0 or x+n>=image.cols or x+n<=0){//do nothing
        }else {
          if(image.at<uchar>(y-j,x+n)>THRESHOLD){imagetemp.at<uchar>(y,x)=0;found=true;};
          if(image.at<uchar>(y+j,x+n)>THRESHOLD){imagetemp.at<uchar>(y,x)=0;found = true;};
        };
      };
      
      j++;    
    /*if(image.at<uchar>(y,x)<200){
      imagetemp.at<uchar>(y,x)=0;
    };*/
    };  
  };
 //cv::namedWindow(window_name1,cv::WINDOW_NORMAL);
 cv::namedWindow(window_name2,cv::WINDOW_NORMAL);
 //cv::imshow(window_name1,image);
 cv::imshow(window_name2,imagetemp);
 cv::waitKey(0);
}