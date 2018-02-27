#ifndef SAC_HPP
#define SAC_HPP
#include<armadillo>
#include <iostream>

struct timeInt{
    double start;
    double end;
};

template <class system, class objective>
class sac {
    system* sys; //from sys use sys->f, sys->proj_func, sys->dfdx
    objective* cost; //from cost use cost->l, cost->dldx, cost->costcalc
    public:
    sac(system *_sys, objective *_cost){
        sys = _sys; cost=_cost;
    };
    
    //algorithm parameters
    double gamma = -5; double delt_init = 0.2; double beta = 0.55;
    double tcalc = 0.0; int kmax = 6; double T = 1.0; 
    arma::vec umax = {20};
    //required functions for calc
    //xforward; rhobackward; MinDisc;
    arma::mat xforward(const arma::vec& u);
    arma::mat rhoback(const arma::mat& xsol);
    
    arma::vec saturation(const arma::vec& u){
        arma::vec usat; usat.zeros(u.n_rows);
        for (int i = 0; i<u.n_rows; i++){
                std::cout<<usat(i)<<"\n";
                if(u(i) > umax(i)) usat(i) = umax(i);
                else if(u(i) < -umax(i)) usat(i) = -umax(i);
                else usat(i) = u(i);
            };
            return usat;
    }
    inline timeInt uInterval(double controlInt[], double simInt[]){
        timeInt interval;
        interval.start= (*simInt); interval.end = *(simInt+1);
        if(controlInt[0]>=simInt[0]) interval.start = controlInt[0];
        else if (controlInt[1]<=simInt[1]) interval.end = controlInt[1];
        return interval;
        };
    
    

};

#endif