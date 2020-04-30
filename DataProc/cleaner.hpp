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
    double scale;
    double offset;
    double mir;
    double refl;
};

class cleaner{
  imagewalls apple("apple.png",BOUND,1.0,1.0);
  imagewalls banana("banana.png",BOUND,1.0,1.0);
  imagewalls house("house.png",BOUND,1.0,1.0);
  imagewalls umbrella("umbrella.png",BOUND,1.0,1.0);
  vector<string> filetypes;
    
  inline int transformcoords(float q, coordtrans Kq){
    return round(Kq.refl+(Kq.mir*Kq.scale*q)+Kq.offset);}
    
      
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
    }
    void crop_data(int s,int f);
    coordtrans transelect(int s, int f, int i);
    //function to select correct image comparison
}

coordtrans cleaner::transelect(int sub, int fnum, int inum){
    coordtrans K;
    K.scale = 2200.;K.offset=0.; K.mir=1.;K.refl=0.;
    if(fnum>0 and fnum <=4){K.scale=1.;}
    if(fnum>4 and fnum <=8){K.offset = 0.5;}
    if(sub<=13){
        if(inum ==1 or inum ==4){K.mir=-1.;K.refl=2200.};
    }
    
}

void cleaner::crop_data(int sub, int fnum){
    //make it work for the simple unmixed image case
    imagewalls walltest(imageName, BOUND,1.0,1.0);
    
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
    //string time,xtemp,ytemp,temp;
    float xt,yt;
    while(getline(myFile,line)){
        vector<string> datastr;
        stringstream ss(line);
        while(ss.good()){
            string substr;
            getline(ss,substr,',');
            datastr.push_back(substr);}
        if(stof(datastr.at(1))<=10.0){
            neighbor nearestpix;
            int x = round(2200.*stof(datastr.at(xind)));
            int y = round(2200.*stof(datastr.at(xind+1)));
            nearestpix = walltest.findnearest(x,y);
            newfile<<line;
            newfile<<","<<nearestpix.dist<<"\n";
        };
    };
    newfile.close();
    
};

#endif