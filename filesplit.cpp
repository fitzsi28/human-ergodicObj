#include "DataProc/cleaner.hpp"
#include<string>
#include<vector>
#include<omp.h>

string DIR = "/home/kt-fitz/data/raw/";
string WRITE = "/home/kt-fitz/data/cpp2/";

int main(){ 
    //omp_set_dynamic(0); // get rid of dynamic stuff
    //omp_set_num_threads(8); // set the number of threads
    cleaner janit(DIR,WRITE);
    //#pragma omp parallel for
    for(int s=2;s<14;s++){
        for(int f=0;f<13;f++){
        //cleaner janit(DIR,WRITE);
        janit.crop_data(s,f);};
    }
    
}
