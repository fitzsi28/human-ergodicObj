#ifndef SAC_HPP
#define SAC_HPP
#include<armadillo>

template <class system, class objective>
class sac {
    system* sys; //from sys use sys->f, sys->proj_func, sys->dfdx
    objective* cost; //from cost use cost->l, cost->dldx, cost->costcalc
    public:
    sac(system *_sys, objective *_cost);
    
    //algorithm parameters
    double gamma = -5; double delt_init = 0.2; double beta = 0.55;
    double tcalc = 0.0; int kmax = 6; double T = 1.0; double umax = 20;
    //required functions for calc
    //xforward; rhobackward;saturation; uInterval; MinDisc;
    arma::mat xforward(const arma::vec& u);
    arma::mat rhoback(const arma::mat& xsol);
    arma::vec saturation(const arma::vec& u);
    inline double * uInterval(double controlInt[], double simInt[]){
        double interval [2] = simInt;
        if(controlInt[0]>=simInt[0]) interval[0] = controlInt[0];
        else if (controlInt[1]<=simInt[1]) interval[1] = controlInt[1];
        return interval;
        };
    
    

};

#endif