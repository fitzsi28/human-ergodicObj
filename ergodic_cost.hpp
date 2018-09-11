#ifndef ERGODICCOST_HPP
#define ERGODICCOST_HPP
#include<armadillo>

template <class system>
class ergodicost {
  system* sys;
  double L1,L2,L1a,L2a, T;
  int X1,X2;//index of relavant dimensions in xvector
  arma::mat hk;
  arma::mat phik;
  arma::mat ckpast;
  arma::mat cktemp;
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
  
  public:
    double Q;
    arma::mat R;
    std::function<double(double,double)> phid;
    ergodicost(double _Q, arma::mat _R,int _K, int _X1,int _X2,std::function<double(double,double)> _phid,double (&_Ln)[2][2],
        double _T,system *_sys){
      Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; T=_T; // initialize with Q, R, sys, phid, and the domain
      X1 = _X1; X2=_X2;
      L1a=_Ln[0][0]; L2a=_Ln[1][0]; L1 =_Ln[0][1]-L1a; L2=_Ln[1][1]-L2a;
      hk.set_size(K,K); hkfunc(); 
      phik.set_size(K,K); phikfunc();
      ckpast.zeros(K,K); cktemp.zeros(K,K);
    };
    
    arma::vec dldx (const arma::vec&x, const arma::vec& u, double ti);
    double calc_cost (const arma::mat& x,const arma::mat& u);
    void ckmemory (const arma::vec&);
};/////////end main class def
template<class system> arma::vec ergodicost<system>::dldx (const arma::vec&x, const arma::vec& u, double ti){
  arma::vec xproj = sys->proj_func(x);
  xproj(X1) = xproj(X1)-L1a; xproj(X2) = xproj(X2)-L2a;
  arma::vec a; a.zeros(xproj.n_rows);
  double LamK;
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),1.5);
      a(X1)+=LamK*2.*(cktemp(k1,k2)-phik(k1,k2))/T*-sin(k1*PI*xproj(X1)/L1)*cos(k2*PI*xproj(X2)/L2)/hk(k1,k2)*k1*PI/L1;
      a(X2)+=LamK*2.*(cktemp(k1,k2)-phik(k1,k2))/T*cos(k1*PI*xproj(X1)/L1)*-sin(k2*PI*xproj(X2)/L2)/hk(k1,k2)*k2*PI/L2;
    };
  };
return Q*a;}

template<class system> double ergodicost<system>::calc_cost (const arma::mat& x,const arma::mat& u){
  double J1 = 0.;
  double LamK;
  cktemp = ckpast+ckfunc(x);
  //arma::mat ck = ckpast+ckfunc(x);
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),1.5); //cktemp (k1,k2)= ckpast(k1,k2)+ck(k1,k2);
      J1+=arma::as_scalar(LamK*pow((cktemp(k1,k2)-phik(k1,k2)),2.));
    };
  };J1 = Q*J1; 
  for (int i = 0; i<x.n_cols; i++){
    J1+=arma::as_scalar(u.col(i).t()*R*u.col(i));
  };
return J1;}

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

template<class system> void ergodicost<system>::phikfunc(){//subtract for 0 to L1
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
  arma::vec xproj;
  arma::mat ck = arma::zeros<arma::mat>(K,K); 
  for(int m=0;m<K;m++){ 
    for(int n=0;n<K;n++){
      for(int j=0; j<x.n_cols;j++){
        xproj = sys->proj_func(x.col(j)); xproj(X1) = xproj(X1)-L1a; xproj(X2) = xproj(X2)-L2a;
        ck(m,n)+=cos(m*PI*xproj(X1)/L1)*cos(n*PI*xproj(X2)/L2)/hk(m,n);
      };
      ck(m,n)=ck(m,n)/(double)x.n_cols;
    };
  };
return ck;}
template<class system> void ergodicost<system>::ckmemory(const arma::vec& x){
  ckpast=ckpast+ckfunc(x);//cout<<ckfunc(x)<<"\n";
}
#endif