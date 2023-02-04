#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
inline double Franke_1 (Eigen::Vector2d & X){
    return exp(-0.25*(pow(9.0*X(0)-2.0,2)+pow(9.0*X(1)-2.0,2)));
}
inline double Franke_2 (Eigen::Vector2d & X){
    return exp((-1.0/49.0)*pow(9.0*X(0)+1.0,2)-0.1*(9.0*X(1)+1.0));
}
inline double Franke_3 (Eigen::Vector2d & X){
    return exp(-0.25*(pow(9.0*X(0)-7.0,2)+pow(9.0*X(1)-3.0,2)));
}
inline double Franke_4 (Eigen::Vector2d & X){
    return exp(-1.0*(pow(9.0*X(0)-4.0,2)+pow(9.0*X(1)-7.0,2)));
}
inline double Franke_u(Eigen::Vector2d X){
    return 0.75*Franke_1(X) + 0.75*Franke_2(X) +0.5*Franke_3(X) - 0.2*Franke_4(X);
}
inline double Franke_dudx(Eigen::Vector2d & X){
    return -3.375*(9.0*X(0)-2.0)*Franke_1(X)-(27.0/98.0)*(9.0*X(0)+1.0)*Franke_2(X)
    -2.25*(9.0*X(0)-7.0)*Franke_3(X) + 3.6*(9.0*X(0)-4)*Franke_4(X);
}
inline double Franke_dudy(Eigen::Vector2d & X){
    return -3.375*(9.0*X(1)-2.0)*Franke_1(X)-0.675*Franke_2(X)
    -2.25*(9.0*X(1)-3.0)*Franke_3(X)+3.6*(9.0*X(1)-7.0)*Franke_4(X);
}
inline double Franke_d2udxdy(Eigen::Vector2d & X){
    return 15.1875*(9.0*X(1)-2.0)*(9.0*X(0)-2.0)*Franke_1(X)+(243.0/980.0)*(9.0*X(0)+1.0)*Franke_2(X)
    +10.125*(9.0*X(0)-7.0)*(9.0*X(1)-3.0)*Franke_3(X)-64.8*(9.0*X(0)-4.0)*(9.0*X(1)-7.0)*Franke_4(X);
}
inline double Franke_d2udx2(Eigen::Vector2d & X){
    return -3.375*(9.0-4.5*pow(9.0*X(0)-2.0,2))*Franke_1(X)
    -(27.0/98.0)*(9.0-(18.0/49.0)*pow(9.0*X(0)+1.0,2))*Franke_2(X)
    -2.25*(9.0-4.5*pow(9.0*X(0)-7.0,2))*Franke_3(X)
    +3.6*(9.0-18.0*pow(9.0*X(0)-4.0,2))*Franke_4(X);
}
inline double Franke_d2udy2(Eigen::Vector2d & X){
    return -3.375*(9.0-4.5*pow((9.0*X(1)-2.0),2))*Franke_1(X)
    +(243.0/400.0)*Franke_2(X) -2.25*(9.0-4.5*pow(9.0*X(1)-3.0,2))*Franke_3(X)
    +3.6*(9.0-18.0*pow(9.0*X(1)-7.0,2))*Franke_4(X);
}
inline double Franke_c(Eigen::Vector2d X){
    return -0.0;
}
inline double Franke_g(Eigen::Vector2d X){
    return Franke_u(X);
}
inline double Franke_p(Eigen::Vector2d X){
    return Franke_u(X);
}
inline Eigen::Vector2d Franke_b(Eigen::Vector2d X){
    Eigen::Vector2d Output(2);
    return Output* 0.0;
}
inline Eigen::Matrix2d Franke_sigma(Eigen::Vector2d X){
    return Eigen::Matrix2d::Identity(2,2) * 1.41421356237;
}
inline double Franke_f(Eigen::Vector2d X){
    return -Franke_d2udx2(X)-Franke_d2udy2(X);
}
inline double Franke_Varphi(Eigen::Vector2d X,Eigen::Vector2d normal){
    return 0.0;
}
inline double Franke_Psi(Eigen::Vector2d X, Eigen::Vector2d normal){
    return 0.0;
}
inline Eigen::Vector2d Franke_F(Eigen::Vector2d X){
    Eigen::Vector2d F(2);
    F(0) = Franke_dudx(X);
    F(1) = Franke_dudy(X);
    F = -Franke_sigma(X).transpose()*F;
    return F;
}
inline bool Stopping_mix(Eigen::Vector2d X){
    //if(fabs(X(0) + 1.0) < 1E-8) return false;
    return true;
}
inline double Franke_RBF(Eigen::Vector2d x , Eigen::Vector2d xj, double c2){
    return 1/(pow((x-xj).norm(),2) + c2);
}
