#ifndef DKLCOST_HPP
#define DKLCOST_HPP
#include<armadillo>
#include<math.h>
//#include<omp.h>

template <class system>
class dklcost {
  system* sys;
  double L1,L2,T;
  int X1,X2;//index of relavant dimensions in xvector
  arma::uvec X_DKL;
  int K;//number of samples in the search domain
  arma::mat sigma; 
      
  public:
    arma::mat xpast,xfuture;
    arma::mat domainsamps;
    arma::vec qs_i,ps_i;
    double Q;
    arma::mat R;
    int T_index,t_now=0;
    std::function<double(const arma::vec&)> phid;
    dklcost(double _Q, arma::mat _R,int _K,arma::mat _sigma, int _X1,int _X2,std::function<double(const arma::vec&)> _phid,double _L1,double _L2,
        double _T,system *_sys){
      Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; T=_T; // initialize with Q, R, sys, phid, and the domain
      X1 = _X1; X2=_X2; L1 = _L1; L2 = _L2; X_DKL<<X1<<X2; sigma=_sigma;
      T_index = T/sys->dt;
      xfuture=arma::zeros<arma::mat>(sys->Xcurr.n_rows,T_index);
      xpast.set_size(sys->Xcurr.n_rows,300/sys->dt);//initialize xpast to hold up to fiveminutes of data
      domainsamps.set_size(2,K);
      qs_i.zeros(K);ps_i.set_size(K);
      resample(); //initiliaze the uniform samples over the domain and the probability of each sample
    };
    double l (const arma::vec& x,const arma::vec& u,double ti);
    arma::vec dldx (const arma::vec&x, const arma::vec& u, double ti);
    double calc_cost (const arma::mat& x,const arma::mat& u);
    void xmemory (const arma::vec&);
    void resample ();
    void qs_disc(const arma::mat& x);
};/////////end main class def

/*template<class system> double dklcost<system>::l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      arma::mat Qtemp = arma::zeros<arma::mat>(xproj.n_rows,xproj.n_rows);
      Qtemp(X1,X1)= pow(xproj(X1)/(L1+(0.1*L1)),8);
      Qtemp(X2,X2) = pow(xproj(X2)/(L2+(0.1*L2)),8);
      return arma::as_scalar((xproj.t()*Qtemp*xproj+u.t()*R*u)/2);
      }*/
template<class system> arma::vec dklcost<system>::dldx (const arma::vec&x, const arma::vec& u, double ti){
  arma::vec xproj = sys->proj_func(x);
  arma::vec a; a.zeros(xproj.n_rows);
  arma::mat Qtemp = arma::zeros<arma::mat>(xproj.n_rows,xproj.n_rows);
  Qtemp(X1,X1)= pow(xproj(X1)/(L1+(0.1*L1)),8);
  Qtemp(X2,X2) = pow(xproj(X2)/(L2+(0.1*L2)),8);
  //a=a+5*Qtemp*xproj;
  //xproj(X1) = xproj(X1)+L1; xproj(X2) = xproj(X2)+L2;
  //double LamK, Dx1F,Dx2F;
  //arma::vec x_dkl = {xproj(X1),xproj(X2)};
  for(int n=0;n<ps_i.n_rows;n++){
    arma::vec s_x = domainsamps(n)-x.elem(X_DKL);//x_dkl;
    a.elem(X_DKL)-=arma::as_scalar(ps_i(n)/qs_i(n))*exp(-0.5*arma::as_scalar(s_x.t()*sigma.i()*s_x))*s_x.t()*sigma.i();
  };
  //a(X1) = a(X1)+atemp(0); a(X2)=a(X2)+atemp(1);
return a;}

template<class system> double dklcost<system>::calc_cost (const arma::mat& x,const arma::mat& u){
  
  double J1 = 0.; 
  arma::mat xjoined = arma::join_rows(xpast.cols(0,t_now),x);
  qs_disc(xjoined);
  for(int n=0;n<K;n++){
    J1+=arma::as_scalar(ps_i(n)*log(ps_i(n))-ps_i(n)*log(qs_i(n)));
  };J1 = Q*J1; //cout<<"ergodic cost "<<J1<<" ";
  /*for (int i = 0; i<x.n_cols; i++){
    arma::vec xproj = sys->proj_func(x.col(i));
    J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt); 
  };//cout<<"total cost "<<J1<<"\n";*/
return J1;}

template<class system> void dklcost<system>::qs_disc(const arma::mat& x){
  for(int n=0;n<qs_i.n_rows;n++){
    qs_i(n) = 0.;
    for(int j=0;j<x.n_cols;j++){
      arma::vec sj = domainsamps.col(n)-sys->proj_func(x.col(j)).elem(X_DKL);
      qs_i(n)+=sys->dt*exp(-0.5*arma::as_scalar(sj.t()*sigma.i()*sj));
    };
  };
  qs_i=normalise(qs_i,1);//normalise the discrete pdf over the samples
}
      
template<class system> void dklcost<system>::resample(){
  //Choose K samples over the domain [[-L1,L1],[-L2,L2]]
  arma::vec domainsize = {2*L1,2*L2};
  domainsamps=arma::diagmat(domainsize)*arma::randu<arma::mat>(2,K);//generate uniform random samples from 0 to 2*L
  domainsamps.each_col() -= (0.5*domainsize);//shift samples to go from -L to L
  for(int n=0;n<ps_i.n_rows;n++){
    ps_i(n) = phid(domainsamps.col(n));
  };
  ps_i=normalise(ps_i,1);//normalise the discrete pdf over the samples
}

template<class system> void dklcost<system>::xmemory(const arma::vec& x){
  xpast.col(t_now)= x;
  t_now++;
  //resample();
}
#endif