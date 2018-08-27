#ifndef SAC_HPP
#define SAC_HPP
#include<armadillo>
#include <iostream>

struct timeInt{
    double start;
    double end;
};

struct SACaction{
    timeInt tau;
    arma::vec u;
};

struct xupair{
    arma::vec x;
    arma::vec u;
};



template <class system, class objective>
class sac {
    system* sys; //from sys use sys->f, sys->proj_func, sys->dfdx, sys->hx
    objective* cost; //from cost use cost->l, cost->dldx, cost->calc_cost
    
    public:
    
    //algorithm parameters
    double gamma = -5; double delt_init = 0.2; double beta = 0.55;
    double tcalc = 0.0; int kmax = 6; double T = 0.1; 
    int T_index;
    arma::vec umax = {20};
    
    sac(system *_sys, objective *_cost){
        sys = _sys; cost=_cost;
        T_index = T/sys->dt;
        
    };
    
    //main function for calculating a single SAC control vector
    SACaction SAC_calc(const arma::vec& InitCon, const arma::mat& u1);
    //required functions for calc
    //xforward; rhobackward; MinDisc;
    arma::mat xforward(const arma::mat& u);//forward simulation of x
    arma::mat rhoback(const arma::mat& xsol,const arma::mat& u); //backward simulation of the adjoint
        inline arma::vec f(const arma::vec& rho, xupair pair){
            return -cost->dldx(pair.x,pair.u) - sys->dfdx(pair.x,pair.u).t()*rho;
            }//f for rho backwards sim
    //saturation funtion for the control        
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
    //restricting the application interval of ustar to within the next time window
    inline timeInt uInterval(double controlInt[], double simInt[]){
        timeInt interval;
        interval.start= (*simInt); interval.end = *(simInt+1);
        if(controlInt[0]>=simInt[0]) interval.start = controlInt[0];
        else if (controlInt[1]<=simInt[1]) interval.end = controlInt[1];
        return interval;
        };
    
    

};
//main function for calculating a single SAC control vector
template <class system, class objective>
SACaction sac<system,objective>::SAC_calc(const arma::vec& InitCon, const arma::mat& u1){
    SACaction ustar;
    ustar.tau.start = 0; ustar.tau.end = 0; 
    ustar.u = arma::zeros<arma::vec>(size(u1.col(0)));
    arma::uword tautemp;
    arma::mat xsol,rhosol;
    arma::mat usched = arma::zeros<arma::mat>(1,T_index);
    arma::mat Jtau = arma::zeros<arma::mat>(1,T_index);
    double J1init,J1new,dJmin,alphad,lambda;
          
    xsol = xforward(u1);
    rhosol = rhoback(xsol, u1);
    J1init = cost->calc_cost(xsol,u1);
    dJmin = 0;//gamma*J1init
    alphad = gamma*J1init;
    arma::vec Lam;
    arma::vec dJdlam;
    for(int i = 0; i<T_index;i++){
        Lam = sys->hx(xsol.col(i)).t()*rhosol.col(i)*rhosol.col(i).t()*sys->hx(xsol.col(i));
        usched.col(i) = (Lam +cost->R).i()*(Lam*u1.col(i) + sys->hx(xsol.col(i)).t()*rhosol.col(i)*alphad);
        dJdlam = rhosol.col(i).t()*(sys->f(xsol.col(i),usched.col(i))-sys->f(xsol.col(i),u1.col(i)));
        Jtau.col(i) =arma::norm(usched.col(i))+dJdlam+pow(i*sys->dt,beta);
        //if (Jtau(i)
    }
    tautemp = Jtau.index_min();
    cout<<sys->tcurr+(sys->dt*(double)tautemp);
    //ustar.u=usched.col(0);
    ustar.u=usched.col(tautemp);
    int k = 0; J1new = 10*J1init;
    while(J1new-J1init>dJmin && k<= kmax){
        lambda = delt_init*pow(beta,k);
        ustar.tau.start = (double)tautemp*sys->dt -(lambda/2);
        ustar.tau.end = (double)tautemp*sys->dt+(lambda/2);
        k++;
    }
    //Initialize k=0 and J1new = A large number
    //While loop for finding application duration
        //lambda = omega^k delta t init
        //
    
    return ustar;
    }
//forward simulation of x
template <class system, class objective>
arma::mat sac<system,objective>::xforward(const arma::mat& u){
    arma::mat xsol = arma::zeros<arma::mat>(4,T_index);
    arma::vec x0 = sys->Xcurr;
   for(int i = 0; i<T_index;i++){
       xsol.col(i)=x0;
       x0 = RK4_step<system,const arma::vec&>(sys,x0,u.col(i),sys->dt);
   }
    
return xsol;
}

//backward simulaiton of the adjoint
template <class system, class objective>
arma::mat sac<system,objective>::rhoback(const arma::mat& xsol,const arma::mat& u){
    arma::mat rhosol = arma::zeros<arma::mat>(4,T_index);
    arma::vec rho0 = sys->Xcurr;
    xupair current;
    rho0.zeros();
   for(int i = T_index-1; i>=0;i--){
       rhosol.col(i)=rho0;
       current.x =xsol.col(i);
       current.u = u.col(i);
       rho0 = RK4_step<sac,xupair>(this,rho0,current,sys->dt);
   }
    
return rhosol;
}


#endif