#ifndef ERRORCOST_HPP
#define ERRORCOST_HPP
#include<armadillo>

template <class system>
class errorcost {

    public:
        arma::mat Q;
        arma::mat R;
        arma::vec xd;
        system& sys;
        errorcost(arma::mat _Q, arma::mat _R, arma::vec _xd, system& _sys){
            Q=_Q; R=_R; xd=_xd; sys=_sys;// initialize with Q, R, sys, xd
            };
        inline double l (const arma::vec& x,const arma::vec& u){
            arma::vec xproj = sys.proj_func(x);
            return arma::as_scalar(((xproj.t()-xd.t())*Q*(xproj-xd)+u.t()*R*u)/2);
        };
        inline arma::vec dldx (const arma::vec& x,const arma::vec& u){
            arma::vec xproj = sys.proj_func(x);
            return Q*(xproj-xd);
        };
        double calc_cost(const arma::vec& x,const arma::vec& u){
            arma::vec xproj;
            double J1 = 0.0;
            for (int i = 0; i<x.n_cols; i++){
                xproj = sys.proj_func(x.col(i));
                J1+=l(xproj,u.col(i));
            }
            return J1;
        }
    
    };

#endif