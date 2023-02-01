#ifndef INTERFAZ_CLASS
#define INTERFAZ_CLASS
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "Nudo_Interior.hpp"
#include "Nudo_Perimetral.hpp"
#include "../BVP/BVP.hpp"
#define AVERAGED
class Interfaz{
    public:
    std::vector<int> posicion_centro;
    bool es_perimeter;
    Eigen::Vector2d centro;
    double radio,c2;
    int ubicacion, interior;
    std::vector<double> parametros_dominio;
    //Nudos sobre la interfaz.
    std::vector<Nudo> nudos_circunferencia;
    //Nudos de otras interfaces que están dentro de la interfaz.
    std::vector<Nudo> nudos_interior;
    //Matriz de interpolacion RBF y su inversa
    Eigen::MatrixXd Psi, iPsi;
    int Inicia_Cuadrado(double RADIO, BVP problema, std::vector<double> PARAMETROS_DOMINIO, 
    double theta_0, std::map<std::string,double> c2, Eigen::Vector2d CENTRO, int indice_posicion_x,
    int indice_posicion_y, int indice_maximo_x,  int indice_maximo_y,  int numero_nudos_total, 
    int primer_nudo);
    double Calcula_Psi(BVP problema, double C2);
    void Update_G(double Y, Eigen::Vector2d X, BVP problema, std::map<int,double> &G);
};
#endif