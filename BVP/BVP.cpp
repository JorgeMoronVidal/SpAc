#include"BVP.hpp"
BVP::BVP(void){
    Reset();
}
void BVP::Reset(void){
    Coeficiente_Escalar U,G,C,F;
    Coeficiente_EscalarN PSI, VARPHI;
    Coeficiente_Vectorial B, GRADIENTE_U;
    Coeficiente_Matricial SIGMA;
    u = U;
    g = G; 
    c = C; 
    f = F;
    psi = PSI;
    varphi = VARPHI;
    b = B;
    gradiente_u = GRADIENTE_U;
    sigma = SIGMA;
    RBF = Inverse_Multiquadric;
    fcardinal_Fourier = Cardinal_Fourier;
}
inline double Inverse_Multiquadric(double theta, double theta_j, double c2){
    return 1/sqrt(pow(theta-theta_j,2)+c2);
}
double Cardinal_Fourier(double theta, std::vector<double> theta_j){
    double aux = 0.0;
    int n = theta_j.size()/2;
    for(int i = 0; i < theta_j.size(); i++) aux += cos((n+i)*(theta-theta_j[i]));
    return aux/(2*n);
}