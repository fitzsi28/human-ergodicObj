#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv

#include"Virtual_Fixt/imagewalls.hpp"
using namespace std;
const int SEARCHRAD = 100;

int main()
{ const char* window_name1 = "Original";
  const char* window_name2 = "WithinBounds";
  //cv::Mat image;
  string imageName("apple.png");
  imagewalls walltest(imageName, SEARCHRAD,1.0,1.0);
  cv::Mat imagetemp = cv::imread(imageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
  random_device rd; mt19937 eng(rd());
  uniform_int_distribution<> distr(0,2200-1);
  for(int i=0;i<=10000;i++){
    int x = distr(eng);
    int y = distr(eng);
    neighbor nearestpix;
    nearestpix = walltest.findnearest(x,y);
    arma::vec forcetest = walltest.wallforce(x,y);
    if(nearestpix.dist>200.){}//imagetemp.at<uchar>(y,x)=0;}
    else{
      if(nearestpix.dist>25.){imagetemp.at<uchar>(y,x)=0;}
      };
  };
 //cv::namedWindow(window_name1,cv::WINDOW_NORMAL);
 cv::namedWindow(window_name2,cv::WINDOW_NORMAL);
 //cv::imshow(window_name1,image);
 cv::imshow(window_name2,imagetemp);
 cv::waitKey(0);
}