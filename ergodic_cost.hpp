#ifndef ERGODICCOST_HPP
#define ERGODICCOST_HPP
#include<armadillo>

template <class system>
class ergodicost {
        system* sys;
        double L1,L2,L1a,L2a;
        arma::mat hk;
        arma::mat phik;
        arma::mat ck;
        int K;
        inline double trapint(std::function<double(double,double)> f,int N1,int N2,double d1,double d2){
            double total = 0.;
            for(int j=0; j<N2;j++){
                for(int k=0;k<N1;k++){
                    total+=(1./4.)*d1*d2*(f(k*d1,j*d2)+f(k*d1,(j+1)*d2)+f((k+1)*d1,j*d2)+f((k+1)*d1,(j+1)*d2));
                }
            }
        return total;};
        void hkfunc();
        void phikfunc();
        arma::mat ckfunc(const arma::mat& );
        arma::vec dldx (const arma::vec&x, const arma::vec& u, double ti);
    public:
        double Q;
        arma::mat R;
        std::function<double(double,double)> phid;
        ergodicost(double _Q, arma::mat _R,int _K, std::function<double(double,double)> _phid,double (&_Ln)[2][2], system *_sys){
            Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; // initialize with Q, R, sys, phid, and the domain
            L1a=_Ln[0][0]; L2a=_Ln[1][0]; L1 =_Ln[0][1]-L1a; L2=_Ln[1][1]-L2a;
             hk.set_size(K,K); hkfunc(); 
             phik.set_size(K,K); phikfunc();
             ck.set_size(K,K);
             
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
};/////////end main class def
template<class system> arma::vec ergodicost<system>::dldx (const arma::vec&x, const arma::vec& u, double ti){
        arma::vec xproj = sys->proj_func(x);
        arma::vec a; a.zeros(K,K);
        double LamK;
        for(int k1=0;k1<K;k1++){
            for(int k2=0;k2<K;k2++){
            LamK = pow(1+(pow(k1,2)+pow(k2,2)),1.5);
            a(0)+=LamK;
            };
        };
return a;}

template<class system> void ergodicost<system>::hkfunc(){
    for(int m=0;m<K;m++){
        for(int n=0;n<K;n++){
            int L1ind = 100; int L2ind = 100;
            double d1 = L1ind/L1;
            double d2 = L2ind/L2;
            auto fk = [&](double x1,double x2){return pow(cos(m*PI*x1/L1),2.)*pow(cos(n*PI*x2/L2),2.);};
            hk(m,n)=pow(trapint(fk,L1ind,L2ind,d1,d2),0.5);      
        };
    };
}

template<class system> void ergodicost<system>::phikfunc(){
    for(int m=0;m<K;m++){
        for(int n=0;n<K;n++){
            int L1ind = 100; int L2ind = 100;
            double d1 = L1ind/L1;
            double d2 = L2ind/L2;
            auto Fk = [&](double x1,double x2){return phid(x1,x2)*cos(m*PI*x1/L1)*cos(n*PI*x2/L2)/hk(m,n);};
            phik(m,n)=trapint(Fk,L1ind,L2ind,d1,d2);      
        };
    };
}
template<class system> arma::mat ergodicost<system>::ckfunc(const arma::mat& x){
    arma::mat cktemp = arma::zeros<arma::mat>(K,K);
    for(int m=0;m<K;m++){
        for(int n=0;n<K;n++){
            for(int j=0; j<x.n_cols;j++){
                cktemp(m,n)+=cos(m*PI*x(0,j)/L1)*cos(n*PI*x(1,j)/L2)/hk(m,n);
            };
            cktemp(m,n)=cktemp(m,n)/x.n_cols;
        };
    };
return cktemp;}
#endif