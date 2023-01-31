#ifndef BVP_CLASS
#define BVP_CLASS
#include <iostream>
#include "Coeficiente_Escalar.hpp"
#include "Coeficiente_Escalar_Numerico.hpp"
#include "Coeficiente_EscalarN.hpp"
#include "Coeficiente_EscalarN_Numerico.hpp"
#include "Coeficiente_Matricial.hpp"
#include "Coeficiente_Matricial_Numerico.hpp"
#include "Coeficiente_Vectorial.hpp"
#include "Coeficiente_Vectorial_Numerico.hpp"
#include "Frontera.hpp"
typedef double (*pRBF)(double, double, double);
typedef double (*pCardinalFourier)(double, std::vector<double>);
inline double Inverse_Multiquadric(double x, double theta_j, double c2);
double Cardinal_Fourier(double theta, std::vector<double> theta_j);
class BVP{
    public:
    Coeficiente_Escalar u, g, c, f;
    Coeficiente_EscalarN psi, varphi;
    Coeficiente_Vectorial b, gradiente_u,F,mu;
    Coeficiente_Matricial sigma;
    Frontera frontera_dominio, frontera_subdominio;
    pRBF RBF;
    pCardinalFourier fcardinal_Fourier;
    BVP(void);
    void Reset(void);
};
#endif