#include "DataProc/cleaner.hpp"
#include<string>
#include<vector>
#include<omp.h>

string DIR = "/home/kt-fitz/data/raw/";
string WRITE = "/home/kt-fitz/data/cpp/";

int main(){
    vector<string> filetypes;
    filetypes.push_back("_c_set01_mda.csv");
    filetypes.push_back("_apple_v_set02_mda.csv");
    filetypes.push_back("_apple_h_set02_mda.csv");
    filetypes.push_back("_apple_p_set03_mda.csv");
    filetypes.push_back("_banana_v_set02_mda.csv");
    filetypes.push_back("_banana_h_set02_mda.csv");
    filetypes.push_back("_banana_p_set03_mda.csv");
    filetypes.push_back("_house_v_set02_mda.csv");
    filetypes.push_back("_house_h_set02_mda.csv");
    filetypes.push_back("_house_p_set03_mda.csv");
    filetypes.push_back("_umbrella_v_set02_mda.csv");
    filetypes.push_back("_umbrella_h_set02_mda.csv");
    filetypes.push_back("_umbrella_p_set03_mda.csv");
    
    //#pragma omp parallel for
    char buffer[17]; sprintf(buffer,"s%02d",1);
    cleaner janit(DIR,WRITE);
    janit.crop_data(sub,fnum);//buffer+filetypes.at(2));
    
}

coordtrans (int sub, int fnum){
    scale =2200.
    offset =0.
}