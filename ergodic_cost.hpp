#ifndef ERGODICCOST_HPP
#define ERGODICCOST_HPP
#include<armadillo>

template <class system>
class ergodicost {
        system* sys;
        double L1,L2,L1a,L2a;
        arma::mat hk;
        int K;
        inline double trapint(std::function<double(double,double)> f,int N1,int N2,double d1,double d2){
            double total = 0.;
            for(int j=0; j<N2;j++){
                for(int k=0;j<N1;k++){
                        total+=(1/4)*d1*d2*(f(k*d1,j*d2)+f(k*d1,(j+1)*d2)+f((k+1)*d1,j*d2)+f((k+1)*d1,(j+1)*d2));
                }
            }
        return total;};
        void hkfunc();
    public:
        double Q;
        arma::mat R;
        std::function<arma::vec(double[])> phid;
        ergodicost(double _Q, arma::mat _R,int _K, std::function<arma::vec(double[])> _phid,double (&_Ln)[2][2], system *_sys){
            Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; // initialize with Q, R, sys, phid, and the domain
            L1a=_Ln[0][0]; L2a=_Ln[1][0]; L1 =_Ln[0][1]-L1a; L2=_Ln[1][0]-L2a;
             hk.set_size(K,K); hkfunc();
            };
        
            double phikfunc(int i,int j){
            double phik=0.;
            
            
        return phik;
        };
        /*    
        inline double l (const arma::vec& x,const arma::vec& u,double ti){
            arma::vec xproj = sys->proj_func(x);
            return arma::as_scalar(((xproj.t()-xd(ti).t())*Q*(xproj-xd(ti))+u.t()*R*u)/2);
        }
        inline arma::vec dldx (const arma::vec& x,const arma::vec& u,double ti){
            arma::vec xproj = sys->proj_func(x);
            return Q*(xproj-xd(ti));
        }
        double calc_cost(const arma::mat& x,const arma::mat& u){
            arma::vec xproj;
            double J1 = 0.0;
            for (int i = 0; i<x.n_cols; i++){
                xproj = sys->proj_func(x.col(i));
                J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt);
            }
            return J1;
        }
    */
    };
template<class system> void ergodicost<system>::hkfunc(){
            for(int m=0;m<K;m++){
                for(int n=0;n<K;n++){
                    hk(m,n)=0.0;
                    int L1ind = 100; int L2ind = 100;
                    double d1 = L1ind/L1;
                    double d2 = L2ind/L2;
                    hk(m,n)=trapint([&](double x1,double x2){pow(cos(m*PI*x1/L1),2)*pow(cos(n*PI*x2/L2),2);},L1ind,L2ind,d1,d2);
                 };
            };
        }

#endif