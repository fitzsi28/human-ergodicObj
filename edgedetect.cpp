#include <iostream>
#include<string>
#include <fstream>
#include<math.h>
//#include<armadillo>
#include <opencv2/opencv.hpp>//namespace cv
using namespace std;
const int SEARCHRAD = 50;

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
 /*for(int j=1;j<4;j++){
      for(int m=-j;m<=j;m++){
      cout<<"("<<m<<","<<-j<<")";
      cout<<"("<<m<<","<<j<<")";
      };
     for(int n=-j+1;n<j;n++){
      cout<<"("<<-j<<","<<n<<")";
      cout<<"("<<j<<","<<n<<")";
      };
    cout<<"\n";
    };*/
  random_device rd; mt19937 eng(rd());
  uniform_int_distribution<> distr(SEARCHRAD,2200-1);
  for(int i=0;i<=100000;i++){
    int x = distr(eng);
    int y = distr(eng);
    //int j = 1;
    bool found = false;
    for(int j=1;j<=SEARCHRAD;j++){//while(found ==false or j>10)
      for(int m=-j;m<=j;m++){
          //cout<<"("<<y+m<<", "<<x-j<<"): ";
          //cout<<(double)image.at<uchar>(y+m,x-j);
      if(image.at<uchar>(y+m,x-j)>200){imagetemp.at<uchar>(y,x)=0;found=true;};
      if(image.at<uchar>(y+m,x+j)>200){imagetemp.at<uchar>(y,x)=0;found = true;};
      };
      for(int n=-j+1;n<j;n++){
        if(image.at<uchar>(y-j,x+n)>200){imagetemp.at<uchar>(y,x)=0;found=true;};
        if(image.at<uchar>(y+j,x+n)>200){imagetemp.at<uchar>(y,x)=0;found = true;};
      };
      
      //j++;    
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