#ifndef CLEANER_HPP
#define CLEANER_HPP
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

using namespace std;

#include "Virtual_Fixt/imagewalls.hpp"
const int BOUND = 100;

struct coordtrans{
    float scale;
    float offset;
    float mir;
    float refl;
};

class cleaner{
  imagewalls apple{"apple.png",BOUND,1.0,1.0};
  imagewalls banana{"banana.png",BOUND,1.0,1.0};
  imagewalls house{"house.png",BOUND,1.0,1.0};
  imagewalls umbrella{"umbrella.png",BOUND,1.0,1.0};
  vector<string> filetypes;
    
  inline double transformcoords(float q, coordtrans Kq){
    double ans;
    ans=Kq.refl+(Kq.mir*Kq.scale*(q+Kq.offset));
  return ans;}
    
      
  public:
    string oldir,newdir;
    cleaner(string _oldir, string _newdir){
      oldir=_oldir;newdir=_newdir;
      filetypes.push_back("_c_set01_mda.csv");
      filetypes.push_back("_apple_v_set02_mda.csv");
      filetypes.push_back("_banana_v_set02_mda.csv");
      filetypes.push_back("_house_v_set02_mda.csv");
      filetypes.push_back("_umbrella_v_set02_mda.csv");
      filetypes.push_back("_apple_h_set02_mda.csv");
      filetypes.push_back("_banana_h_set02_mda.csv");
      filetypes.push_back("_house_h_set02_mda.csv");
      filetypes.push_back("_umbrella_h_set02_mda.csv");
      filetypes.push_back("_apple_p_set03_mda.csv");
      filetypes.push_back("_banana_p_set03_mda.csv");
      filetypes.push_back("_house_p_set03_mda.csv");
      filetypes.push_back("_umbrella_p_set03_mda.csv");
    };
    void crop_data(int s,int f);
    coordtrans transelect(int s, int f, int i);
    int imageselect(int f);
    //function to select correct image comparison
};

coordtrans cleaner::transelect(int sub, int fnum, int inum){
    coordtrans K;
    K.scale = 2200;K.offset=0.; K.mir=1.;K.refl=0.;
    if(fnum>0 and fnum <=4){K.scale=1;}
    if(fnum>4 and fnum <=8){K.offset = 0.5;}
    if(sub<=13){
        if(inum ==1 or inum ==4){K.mir=-1.;K.refl=2200.;};
    }
  return K;  
}
int cleaner::imageselect(int fnum){
    int i=0;
    if(fnum==1 or fnum ==5 or fnum==9){i=1;}
    else if(fnum==2 or fnum==6 or fnum==10){i=2;}
    else if(fnum==3 or fnum==7 or fnum==11){i=3;}
    else if(fnum==4 or fnum==8 or fnum==12){i=4;}
    return i;    
}
void cleaner::crop_data(int sub, int fnum){
    vector<string> imlabels;
    int imnum = imageselect(fnum);
    coordtrans Kx = transelect(sub,fnum,imnum);
    coordtrans Ky = Kx; Ky.mir=1.;Ky.refl=0.;
    char buffer[17]; sprintf(buffer,"s%02d",sub);
    ifstream imfile(oldir+"image-testing-order-"+buffer+".csv";
    if(imfile.good()==false){return;};
    if(fnum==0){
      string trial;
      while(getline(imfile,trial)){
         stringstream ss(trial);
         while(ss.good()){
           string substr;
           getline(ss,substr,',');
           imlabels.push_back(substr);};
      };
    };
    
    string oldfile = buffer+filetypes.at(fnum);
    string::size_type name_ind = oldfile.rfind('t')+3;
    string filename=newdir+oldfile.substr(0,name_ind)+"-clean.csv";
    ifstream myFile(oldir+oldfile);
    if (myFile.good()==false){return;};
    string colnames,line;
    getline(myFile,colnames);
    int xind = 0;
    string name;
    stringstream cols(colnames);
    int i = 0;
    while(getline(cols,name,',') and xind==0){
        if(name=="field.q0"){xind = i;};
        i++;
    };
    ofstream newfile(filename);
    newfile<<colnames;
    newfile<<",pixdist\n";
    while(getline(myFile,line)){
        vector<string> datastr;
        stringstream ss(line);
        while(ss.good()){
            string substr;
            getline(ss,substr,',');
            datastr.push_back(substr);}
        if(stof(datastr.at(1))<=10.0){
            neighbor nearestpix;
            int x = transformcoords(stof(datastr.at(xind)),Kx);
            int y = transformcoords(stof(datastr.at(xind+1)),Ky);
            switch(imnum){
              case 1: nearestpix = apple.findnearest(x,y);
              case 2: nearestpix = banana.findnearest(x,y);
              case 3: nearestpix = house.findnearest(x,y);
              case 4: nearestpix = umbrella.findnearest(x,y);
            }
            newfile<<line;
            newfile<<","<<nearestpix.dist<<"\n";
        };
    };
    newfile.close();
    
};

#endif