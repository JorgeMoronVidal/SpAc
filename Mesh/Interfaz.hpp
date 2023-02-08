#ifndef INTERFAZ_CLASS
#define INTERFAZ_CLASS
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "Nudo.hpp"
#include "../BVP/BVP.hpp"
#define AVERAGED
class Interfaz{
    public:
    bool es_perimeter;
    std::vector<int> posicion_centro;
    int ubicacion, interior;
    Eigen::Vector2d centro;
    double radio,c2;
    std::vector<double> parametros_dominio;
    //Matriz de interpolacion RBF y su inversa
    Eigen::MatrixXd Psi, iPsi;
    //Nudos sobre la interfaz.
    std::vector<Nudo> nudos_circunferencia;
    //Nudos de otras interfaces que están dentro de la interfaz.
    std::vector<Nudo> nudos_interior;
    int Inicia_Cuadrado(double RADIO, BVP problema, std::vector<double> PARAMETROS_DOMINIO, 
    double theta_0, std::map<std::string,double> c2, Eigen::Vector2d CENTRO, int indice_posicion_x,
    int indice_posicion_y, int indice_maximo_x,  int indice_maximo_y,  int numero_nudos_total, 
    int primer_nudo);
    double Calcula_Psi(BVP problema, double C2);
    void Update_G(double Y, Eigen::Vector2d X, BVP problema, std::map<int,double> &G);
    void Resetea(void);
    void Empaca(std::vector<int> & Interfaz_int, std::vector<double> & Interfaz_double, std::vector<Nudo> &Interfaz_nudos);
    void Desempaca(std::vector<int> & Interfaz_int, std::vector<double> & Interfaz_double, std::vector<Nudo> &Interfaz_nudos);
};
#endif