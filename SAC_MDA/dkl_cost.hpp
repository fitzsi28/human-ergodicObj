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
  int K;//number of samples in the search domain
  inline double trapint(std::function<double(double,double)> f){
    double total = 0.,x0,y0,xf,yf;
    //x0=-L1,y0=-L2,xf=L1,yf=L2;
    x0=0.,y0=0.,xf=2*L1,yf=2*L2;
    int N1=50,N2=50; double d1=(2.*L1/N1);double d2=(2.*L2/N2);
    total = (d1*d2/4)*(f(x0,y0)+f(x0,yf)+f(xf,y0)+f(xf,yf));
    for(int i = 1; i<N1;i++){
      total+=(d1*d2/2)*(f(x0+i*d1,y0)+f(x0+i*d1,yf));
    }
    for(int j = 1; j<N2;j++){
      total+=(d1*d2/2)*(f(x0,y0+j*d2)+f(xf,y0+j*d2));
    }
    for(int k=1; k<N2;k++){
      for(int j=1;j<N1;j++){ 
        total+=d1*d2*f(x0+k*d1,y0+j*d2);
      }
    }
  return total;};
  /*
  inline double eulint(const arma::mat& x,int m, int n){
    arma::vec xproj;
    double total = 0.;
    for(int j=0; j<x.n_cols;j++){
      xproj = sys->proj_func(x.col(j)); xproj(X1) = xproj(X1)+L1; xproj(X2) = xproj(X2)+L2;
      total+=sys->dt*cos(m*PI*xproj(X1)/(2*L1))*cos(n*PI*xproj(X2)/(2*L2))/hk(m,n); 
    };
    return total;};*/
    
  public:
    arma::mat xpast;
    arma::mat domainsamps;
    double Q;
    arma::mat R;
    std::function<double(double,double)> phid;
    dklcost(double _Q, arma::mat _R,int _K, int _X1,int _X2,std::function<double(double,double)> _phid,double _L1,double _L2,
        double _T,system *_sys){
      Q=_Q; R=_R; sys=_sys; K = _K; phid = _phid; T=_T; // initialize with Q, R, sys, phid, and the domain
      X1 = _X1; X2=_X2; L1 = _L1; L2 = _L2;
      //Choose K samples over the domain [[0,2*L1],[0,2*L2]]
      
    };
    double l (const arma::vec& x,const arma::vec& u,double ti);
    arma::vec dldx (const arma::vec&x, const arma::vec& u, double ti);
    double calc_cost (const arma::mat& x,const arma::mat& u);
    void xmemory (const arma::vec&);
};/////////end main class def

template<class system> double dklcost<system>::l (const arma::vec& x,const arma::vec& u,double ti){
      arma::vec xproj = sys->proj_func(x);
      arma::mat Qtemp = arma::zeros<arma::mat>(xproj.n_rows,xproj.n_rows);
      Qtemp(X1,X1)= pow(xproj(X1)/(L1+(0.1*L1)),8);
      Qtemp(X2,X2) = pow(xproj(X2)/(L2+(0.1*L2)),8);
      return arma::as_scalar((xproj.t()*Qtemp*xproj+u.t()*R*u)/2);
      }
template<class system> arma::vec dklcost<system>::dldx (const arma::vec&x, const arma::vec& u, double ti){
  arma::vec xproj = sys->proj_func(x);
  arma::vec a; a.zeros(xproj.n_rows);
  arma::mat Qtemp = arma::zeros<arma::mat>(xproj.n_rows,xproj.n_rows);
  Qtemp(X1,X1)= pow(xproj(X1)/(L1+(0.1*L1)),8);
  Qtemp(X2,X2) = pow(xproj(X2)/(L2+(0.1*L2)),8);
  a=a+5*Qtemp*xproj;
  xproj(X1) = xproj(X1)+L1; xproj(X2) = xproj(X2)+L2;
  double LamK, Dx1F,Dx2F;
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),-1.5);
      Dx1F = arma::as_scalar((-k1*PI/(2*L1*hk(k1,k2)))*sin(k1*PI*xproj(X1)/(2*L1))*cos(k2*PI*xproj(X2)/(2*L2)));
      Dx2F = arma::as_scalar((-k2*PI/(2*L2*hk(k1,k2)))*cos(k1*PI*xproj(X1)/(2*L1))*sin(k2*PI*xproj(X2)/(2*L2)));
      a(X1)+=Q*LamK*2.*(cktemp(k1,k2)-phik(k1,k2))/T*Dx1F;
      a(X2)+=Q*LamK*2.*(cktemp(k1,k2)-phik(k1,k2))/T*Dx2F;

    };
  };
return a;}

template<class system> double dklcost<system>::calc_cost (const arma::mat& x,const arma::mat& u){
  
  double J1 = 0.; double LamK;
  cktemp = (ckpast*sys->tcurr/(sys->tcurr+T))+ckfunc(x);
  for(int k1=0;k1<K;k1++){
    for(int k2=0;k2<K;k2++){
      LamK = pow(1+(pow(k1,2)+pow(k2,2)),-1.5); 
      J1+=arma::as_scalar(LamK*pow((cktemp(k1,k2)-phik(k1,k2)),2.));
    };
  };J1 = Q*J1; //cout<<"ergodic cost "<<J1<<" ";
  for (int i = 0; i<x.n_cols; i++){
    arma::vec xproj = sys->proj_func(x.col(i));
    J1+=l(xproj,u.col(i),sys->tcurr+(double)i*sys->dt); 
  };//cout<<"total cost "<<J1<<"\n";
return J1;}

template<class system> double dklcost<system>::gpast(const arma::mat& s){//g(s) defined for 0 to tcurr
  double total = 0.;
  int t_index = round(sys->tcurr/sys->dt);
  for(int i = 0;i<t_index;i++{
    arma::vec si = s-sys->proj_func(xpast.col(i));
    total+=sys->dt*exp(-0.5*arma::as_scalar(si.t()*sigma.i()*si));
  };
}

template<class system> void ergodicost<system>::phikfunc(){//integrate for 0 to L1
  cout<<trapint(phid)<<"\n";
  #pragma omp parallel for 
  for(int n=0;n<K;n++){
    //#pragma omp parallel for
    for(int m=0;m<K;m++){
      int L1ind = 100; int L2ind = 100;
      double d1 = 2*L1/L1ind;
      double d2 = 2*L2/L2ind;
      auto Fk = [&](double x1,double x2){
        return phid(x1,x2)*cos(m*PI*(x1)/(2*L1))*cos(n*PI*(x2)/(2*L2))/hk(m,n);};
      phik(m,n)=trapint(Fk);      
    };
  };
}
template<class system> arma::mat ergodicost<system>::ckfunc(const arma::mat& x){
  arma::vec xproj;
  arma::mat ck = arma::zeros<arma::mat>(K,K); 
  #pragma omp parallel for 
  for(int m=0;m<K;m++){
    //#pragma omp parallel for
    for(int n=0;n<K;n++){
      ck(m,n) = eulint(x,m,n);
      ck(m,n)=ck(m,n)/(sys->tcurr+x.n_cols*sys->dt);
    };
  };
return ck;}
template<class system> void dklcost<system>::xmemory(const arma::vec& x){
  ckpast=(ckpast*sys->tcurr/(sys->tcurr+sys->dt))+ckfunc(x);//cout<<ckfunc(x)<<"\n";
}
#endif