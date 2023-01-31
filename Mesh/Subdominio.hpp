#ifndef SUBDOMINIO_CLASS
#define SUBDOMINIO_CLASS
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <eigen3/Eigen/Core>
#include "../BVP/BVP.hpp"
class Subdominio{
    public:
    std::vector<int> posicion;
    bool es_perimeter;
    std::map<std::string,TablaDeValores> datos;
    Eigen::Vector2d centro;
    double radio;
    int iteracion;
    int ubicacion;
    Subdominio(void){
        //To be coded
    };
};
#endif