#include <iostream>
#include<math.h>
using namespace std;
const double PI = 3.1415926535987;
float AngWrap (float th);


int main()
{
    const float m=0.1, B = 0.01, g = 9.81, h=2.0;//System Parameters
    float th;
    th = -6.1; 
    cout<<"Hello World!\n";
    cout<<AngWrap(th);cout<<"\n";
}

float AngWrap (float th){
    float thwrap;
    thwrap = fmod(th+PI, 2*PI);
    if (thwrap < 0.0) thwrap = thwrap + 2*PI;
    thwrap = thwrap - PI;
    return thwrap;
}