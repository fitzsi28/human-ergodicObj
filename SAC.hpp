#ifndef SAC_HPP
#define SAC_HPP
#include<armadillo>

template <class system, class objective>
class sac {

    sac(system& sys, objective& cost);
    //from sys use f, proj_func, dfdx
    //from cost use l, dldx, costcalc 
    //algorithm parameters
    //gamma; delt_init;beta;tcalc;kmax;
    //required functions for calc
    //xforward; rhobackward; costcalc;saturation; uInterval; MinDisc;

};

#endif