#ifndef SAC_HPP
#define SAC_HPP
#include<armadillo>
#include <iostream>

struct timeInt{
    double start;
    double end;
};

struct xupair{
    arma::vec x;
    arma::vec u;
};

/*
template<class system, class objective>
class rhoclass{
    system* sys;
    objective* cost;
    //arma::vec rho;
    public:
    rhoclass(system *_sys, objective *_cost){
    sys = _sys; cost=_cost;  
    }
    
    
};*/

template <class system, class objective>
class sac {
    system* sys; //from sys use sys->f, sys->proj_func, sys->dfdx
    objective* cost; //from cost use cost->l, cost->dldx, cost->costcalc
    //rhoclass<system,objective> rhodot;//(sys,cost);
    
    public:
    
    //algorithm parameters
    double gamma = -5; double delt_init = 0.2; double beta = 0.55;
    double tcalc = 0.0; int kmax = 6; double T = 0.05; 
    int T_index;
    arma::vec umax = {20};
    
    sac(system *_sys, objective *_cost){
        sys = _sys; cost=_cost;
        T_index = T/sys->dt;
        
    };
    
    
    //required functions for calc
    //xforward; rhobackward; MinDisc;
    arma::mat xforward(const arma::vec& u);
    arma::mat rhoback(const arma::mat& xsol,const arma::vec& u);
        inline arma::vec f(const arma::vec& rho, xupair pair){
            return -cost->dldx(pair.x,pair.u) - sys->dfdx(pair.x,pair.u).t()*rho;
            }//f for rho backwards sim
            
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

template <class system, class objective>
arma::mat sac<system,objective>::xforward(const arma::vec& u){
    arma::mat xsol = arma::zeros<arma::mat>(4,T_index);
    arma::vec x0 = sys->Xcurr;
   for(int i = 0; i<T_index;i++){
       xsol.col(i)=x0;
       x0 = RK4_step<system,const arma::vec&>(sys,x0,u,sys->dt);
   }
    
return xsol;
}

template <class system, class objective>
arma::mat sac<system,objective>::rhoback(const arma::mat& xsol,const arma::vec& u){
    arma::mat rhosol = arma::zeros<arma::mat>(4,T_index);
    arma::vec rho0 = sys->Xcurr;
    xupair current;
    rho0.zeros();
   for(int i = T_index-1; i>0;i--){
       rhosol.col(i)=rho0;
       current.x =xsol.col(i);
       current.u = u;
       rho0 = RK4_step<sac,xupair>(this,rho0,current,sys->dt);
   }
    
return rhosol;
}


#endif