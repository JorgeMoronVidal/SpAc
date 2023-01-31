#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
#define k_1 1.0
#define k_2 0.01 
#define k_3 0.02 
#define k_4 0.12
#define k_5 0.05
#define k_6 0.05 
#define k_7 -0.12
#define scaling_factor 0.333
inline double Monegros_u(Eigen::Vector2d X){
    return 3.0 + scaling_factor*(sin(sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1))) + tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1))));
}
inline double Monegros_dudx(Eigen::Vector2d X){
    return scaling_factor*((cos(sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)))*k_2*X(0)/sqrt(k_1 + k_2*X(0)*X(0)+k_3*X(1)*X(1))) -
    (k_4*cos(k_4*X(0)+k_5*X(1))+k_6*cos(k_6*X(0)+k_7*X(1)))*(pow(tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1))),2)-1));
}
inline double Monegros_dudy(Eigen::Vector2d X){
    return scaling_factor*(cos(sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)))*k_1*k_3*X(1)/sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)) -
    (k_5*cos(k_4*X(0)+k_5*X(1))+k_7*cos(k_6*X(0)+k_7*X(1)))*(pow(tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1))),2)-1));
}
inline double Monegros_d2udx2(Eigen::Vector2d X){
    double aux1 = sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)),aux2 = tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1)));
    return scaling_factor*(((k_2/aux1)-(k_2*k_2*X(0)*X(0)/pow(aux1,3)))*cos(aux1)-(X(0)*X(0)*k_2*k_2/(aux1*aux1))*sin(aux1)
    +2.0*pow(k_4*cos(k_4*X(0)+k_5*X(1))+k_6*cos(k_6*X(0)+k_7*X(1)),2)*(pow(aux2,2)-1)*aux2
    +(k_4*k_4*sin(k_4*X(0)+k_5*X(1))+k_6*k_6*sin(k_6*X(0)+k_7*X(1)))*(pow(aux2,2)-1));
}
inline double Monegros_d2udy2(Eigen::Vector2d X){
    double aux1 = sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)),aux2 = tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1)));
    return scaling_factor*(((k_3/aux1)-(k_3*k_3*X(1)*X(1)/pow(aux1,3)))*cos(aux1)-(X(1)*X(1)*k_3*k_3/(aux1*aux1))*sin(aux1)
    +2.0*pow(k_5*cos(k_4*X(0)+k_5*X(1))+k_7*cos(k_6*X(0)+k_7*X(1)),2)*(pow(aux2,2)-1)*aux2
    +(k_5*k_5*sin(k_4*X(0)+k_5*X(1))+k_7*k_7*sin(k_6*X(0)+k_7*X(1)))*(pow(aux2,2)-1));
}
inline double Monegros_d2udxdy(Eigen::Vector2d X){
    double aux1 = sqrt(k_1+k_2*X(0)*X(0)+k_3*X(1)*X(1)),aux2 = tanh(sin(k_4*X(0)+k_5*X(1))+sin(k_6*X(0)+k_7*X(1)));
    return scaling_factor*(-(k_2*k_3*X(0)*X(1)*sin(aux1)/(aux1*aux1)) - (k_2*k_3*X(0)*X(1)/pow(aux1,3))*cos(aux1) +
    2.0*(k_4*cos(k_4*X(0)+k_5*X(1))+k_6*cos(k_6*X(0)+k_7*X(1)))*(k_5*cos(k_4*X(0)+k_5*X(1))+k_7*cos(k_6*X(0)+k_7*X(1)))*
    (aux2*aux2 -1)*aux2 + (k_4*k_5*sin(k_4*X(0)+k_5*X(1)) +  k_6*k_7*sin(k_6*X(0)+k_7*X(1)))*(aux2*aux2-1));
}
inline double Monegros_c(Eigen::Vector2d X){
    return -0.0;
}
inline double Monegros_g(Eigen::Vector2d X){
    return Monegros_u(X);
}
inline double Monegros_p(Eigen::Vector2d X){
    return Monegros_u(X);
}
inline Eigen::Vector2d Monegros_b(Eigen::Vector2d X){
    Eigen::Vector2d Output(2);
    return Output* 0.0;
}
inline Eigen::Matrix2d Monegros_sigma(Eigen::Vector2d X){
    return Eigen::Matrix2d::Identity() * 1.41421356237;
}
inline double Monegros_f(Eigen::Vector2d X){
    return -Monegros_d2udx2(X)-Monegros_d2udy2(X);
}
inline double Monegros_f_NL(Eigen::Vector2d X){
    return -Monegros_d2udx2(X)-Monegros_d2udy2(X)+pow(Monegros_u(X),3);
}
inline double Monegros_Varphi(Eigen::Vector2d X,Eigen::Vector2d normal){
    return 0.0;
}
inline double Monegros_Psi(Eigen::Vector2d X, Eigen::Vector2d normal){
    return 0.0;
}
inline Eigen::Vector2d Monegros_F(Eigen::Vector2d X){
    Eigen::Vector2d F;
    F(0) = Monegros_dudx(X);
    F(1) = Monegros_dudy(X);
    F = -Monegros_sigma(X).transpose()*F;
    return F;
}
inline Eigen::Vector2d Monegros_grad(Eigen::Vector2d X){
    Eigen::Vector2d grad;
    grad(0) = Monegros_dudx(X);
    grad(1) = Monegros_dudy(X);
    return grad;
}
inline bool Stopping_mix(Eigen::Vector2d X){
    //if(fabs(X(0) + 1.0) < 1E-8) return false;
    return true;
}
inline double Monegros_RBF(Eigen::Vector2d x , Eigen::Vector2d xj, double c2){
    return 1/(pow((x-xj).norm(),2) + c2);
}