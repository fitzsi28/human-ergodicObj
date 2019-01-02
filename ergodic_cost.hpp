#ifndef ERGODICCOST_HPP
#define ERGODICCOST_HPP
#include<armadillo>
#include<math.h>

template <class system>
class ergodicost {
  system* sys;
  double L1,L2,T;
  int X1,X2;//index of relavant dimensions in xvector
  arma::mat hk;
  arma::mat phik;
  arma::mat ckpast;
  arma::mat cktemp;
  int K;
  inline double trapint(std::function<double(double,double)> f,int N1,int N2,double d1,double d2,double x0,double y0){
    double total = 0.;double x1; double x2;
    for(int j=0; j<N2;j++){
      for(int k=0;k<N1;k++){ x1 = x0 +k*d1; x2 = y0+j*d1;
        total+=(1./4.)*d1*d2*(f(k*d1,j*d2)+f(k*d1,(j+1)*d2)+f((k+1)*d1,j*d2)+f((k+1)*d1,(j+1)*d2));
      }
    }
  return total;};
  void hkfunc();
  void phikfunc();
  arma::mat ckfunc(const arma::mat& );
  arma::mat Q2 = {
        {0.,0.,0.,0.},
        {0., 20.,0.,0.},
        {0.,0.,40.,0.},
        {0.,0.,0.,1.}};
  public:
    double Q;
    arma::mat R;
    std::function<double(double,double)> phid;
    ergodicost(double _Q, arma::mat _R,int _K, int _X1,int _X2,std::function<double(double,double)> _phid,double _L1,double _L2,
        double _T,system *_sys){
      Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; T=_T; // initialize with Q, R, sys, phid, and the domain
      X1 = _X1; X2=_X2; L1 = _L1; L2 = _L2;
      hk.set_size(K,K); hkfunc(); 
      phik.set_size(K,K); phikfunc();
      ckpast.zeros(K,K); cktemp.zeros(K,K);
    };
    double l (const arma::vec& x,const arma::vec& u,double ti);
    arma::vec dldx (const arma::vec&x, const arma::vec& u, double ti);
    double calc_cost (const arma::mat& x,const arma::mat& u);
    void ckmemory (const arma::vec&);
};/////////end main class def

template<class system> double ergodicost<system>::l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      return arma::as_scalar((xproj.t()*Q2*xproj+u.t()*R*u)/2);
      }
template<class system> arma::vec ergodicost<system>::dldx (const arma::vec&x, const arma::vec& u, double ti){
  arma::vec xproj = sys->proj_func(x);
  arma::vec a; a.zeros(xproj.n_rows);
  a=a+Q2*xproj;
  xproj(X1) = xproj(X1)+L1; xproj(X2) = xproj(X2)+L2;
  double LamK, Dx1F,Dx2F;
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),1.5);
      Dx1F = arma::as_scalar((-k1*PI/(2*L1*hk(k1,k2)))*sin(k1*PI*xproj(X1)/(2*L1))*cos(k2*PI*xproj(X2)/(2*L2)));
      Dx2F = arma::as_scalar((-k2*PI/(2*L2*hk(k1,k2)))*cos(k1*PI*xproj(X1)/(2*L1))*sin(k2*PI*xproj(X2)/(2*L2)));
      a(X1)+=Q*LamK*(cktemp(k1,k2)-phik(k1,k2))/T*Dx1F;
      a(X2)+=Q*LamK*2.*(cktemp(k1,k2)-phik(k1,k2))/T*Dx2F;

    };
  };
return a;}

template<class system> double ergodicost<system>::calc_cost (const arma::mat& x,const arma::mat& u){
  
  double J1 = 0.; double LamK;
  cktemp = ckpast+ckfunc(x);
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),1.5); 
      J1+=arma::as_scalar(LamK*pow((cktemp(k1,k2)-phik(k1,k2)),2.));
    };
  };J1 = Q*J1; 
  for (int i = 0; i<x.n_cols; i++){
    arma::vec xproj = sys->proj_func(x.col(i));
    J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt);
  };
return J1;}

template<class system> void ergodicost<system>::hkfunc(){//integrate 0 to L
  for(int m=0;m<K;m++){
    for(int n=0;n<K;n++){
      int L1ind = 50; int L2ind = 50;
      double d1 = 2*L1/L1ind;
      double d2 = 2*L2/L2ind;
      auto fk = [&](double x1,double x2){
        return pow(cos(m*PI*(x1+L1)/(2*L1)),2.)*pow(cos(n*PI*(x2+L2)/(2*L2)),2.);};
      hk(m,n)=pow(trapint(fk,L1ind,L2ind,d1,d2,0.,0.),0.5);      
    };
  };
}

template<class system> void ergodicost<system>::phikfunc(){//integrate for 0 to L1
  for(int m=0;m<K;m++){
    for(int n=0;n<K;n++){
      int L1ind = 50; int L2ind = 50;
      double d1 = 2*L1/L1ind;
      double d2 = 2*L2/L2ind;
      auto Fk = [&](double x1,double x2){
        return phid(x1,x2)*cos(m*PI*(x1+L1)/(2*L1))*cos(n*PI*(x2+L2)/(2*L2))/hk(m,n);};
      phik(m,n)=trapint(Fk,L1ind,L2ind,d1,d2,0.,0.);      
    };
  };
}
template<class system> arma::mat ergodicost<system>::ckfunc(const arma::mat& x){
  arma::vec xproj;
  arma::mat ck = arma::zeros<arma::mat>(K,K); 
  for(int m=0;m<K;m++){ 
    for(int n=0;n<K;n++){
      for(int j=0; j<x.n_cols;j++){
        xproj = sys->proj_func(x.col(j)); xproj(X1) = xproj(X1)+L1; xproj(X2) = xproj(X2)+L2;
        ck(m,n)+=cos(m*PI*xproj(X1)/(2*L1))*cos(n*PI*xproj(X2)/(2*L2))/hk(m,n);
      };
      ck(m,n)=ck(m,n)/(double)x.n_cols;
    };
  };
return ck;}
template<class system> void ergodicost<system>::ckmemory(const arma::vec& x){
  ckpast=ckpast+ckfunc(x);//cout<<ckfunc(x)<<"\n";
}
#endif