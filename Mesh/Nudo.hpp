#ifndef NUDO_CLASS
#define NUDO_CLASS
#include <iostream>
#include <eigen3/Eigen/Core>
#include <vector>
#include <map>
inline double Atan0(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux >= 0) return aux;
    return aux + 2.0*M_PI;
}
inline double Atan90(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux <= 0.5*M_PI) return 2.0*M_PI + aux;
    return aux;
    
}
inline double Atan180(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    return atan2(point(1),point(0));
}
inline double Atan270(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux <= -0.5*M_PI) return 2*M_PI +  aux;
    return aux;
}
class Nudo{
    public:
    //True si el nudo está sobre la frontera, false en cualquier otro caso
    bool frontera;
    //Posiciones cartesiana, angular y dentro de la matriz.
    Eigen::Vector2d posicion_cartesiana;
    double posicion_angular;
    //El indice se identifica con la fila de la matriz[vector] que contiene sus elementos de G[B]
    int indice_global;
    //Indice local dentro de su interfaz
    int indice_local;
    //Subdominios a los que pertenece
    std::vector<std::vector<int>> indice_subdominios;
    //Interfaz a la que pertenece
    std::vector<int> indice_interfaz;
    //Coeficiente B[i], valor de la solucion sobre el nodo, de la varianza, coeficiente de pearson entre phi y xi, 
    //error estadístico cometido y bias.
    double B,solucion, error, varianza, coeficiente_Pearson, error_estadistico, bias;
    //Componentes de G (Formato COO)
    std::map<int,double> G;
    Nudo(void){
        
    }
    virtual void Inicializa(void){
        std::cout << "Nudo. To code \n" << std::endl;
    };
};
#endif