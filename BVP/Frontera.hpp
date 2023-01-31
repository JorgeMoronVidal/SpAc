#ifndef FRONTERA_CLASS
#define FRONTERA_CLASS
#include<iostream>
#include<vector>
#include <eigen3/Eigen/Core>
typedef double (*pfdistancia)(std::vector<double>, 
                             Eigen::Vector2d, 
                             Eigen::Vector2d&,
                             Eigen::Vector2d&);
typedef bool (*pftipo)(Eigen::Vector2d);
//Funcion distancia defecto que hace todas las magnitudes implicadas 0
inline double Distancia_Defecto(std::vector<double> parameters,Eigen::Vector2d position, Eigen::Vector2d & normal, Eigen::Vector2d & nProjection){
        normal = 0.0*position;
        nProjection = normal;
        return 0.0;
}
inline bool Tipo_Defecto(Eigen::Vector2d position){
        return true;
}
class Frontera{
    public:
    pfdistancia Distancia;
    pftipo Absorbente, Neumann;
    std::vector<double> parametros;
    //Constructor por defecto
    Frontera(void);
    //Constructor dada una función 
    Frontera(std::vector<double> parameters, pfdistancia DISTANCIA, pftipo ABSORBENTE, pftipo NEUMANN);
    //Inicializa el objeto dada una función
    virtual void Inicializa(std::vector<double> parameters, pfdistancia DISTANCIA, pftipo ABSORBENTE, pftipo NEUMANN);
    //Devuelve el valor de la función distancia en punto determinado
    virtual double Evalua_Distancia(Eigen::Vector2d punto, Eigen::Vector2d & punto_mas_cercano, Eigen::Vector2d & normal);
    //Verdadero si las condiciones de contorno en el punto evaluado son absorbentes, falso si reflejantes
    bool Evalua_Absorbente(Eigen::Vector2d punto);
    //Verdadero si las condiciones de contorno en el punto evaluado son Neumann, falso si Robin
    bool Evalua_Neumann(Eigen::Vector2d punto);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    virtual bool Es_Numerico(void);
};
#endif