#ifndef COEFICIENTE_ESCALAR
#define COEFICIENTE_ESCALAR
#include<iostream>
#include <eigen3/Eigen/Core>
typedef double (*pfescalar)(Eigen::Vector2d);
//Funcione escalar por defecto que ante cualquier input devuelve siempre 0
inline double Escalar_Defecto(Eigen::Vector2d position){
        return 0.0;
}
class Coeficiente_Escalar{
    public:
    pfescalar funcion;
    //Constructor por defecto
    Coeficiente_Escalar(void);
    //Constructor dada una función 
    Coeficiente_Escalar(pfescalar FUNCION);
    //Inicializa el objeto dada una función
    void Inicializa(pfescalar FUNCION);
    //Devuelve el valor de la función cargada en punto e instante determinados
    double Evalua(Eigen::Vector2d punto);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    bool Es_Numerico(void);
};
#endif