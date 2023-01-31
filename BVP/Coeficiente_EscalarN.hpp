#ifndef COEFICIENTE_ESCALARN
#define COEFICIENTE_ESCALARN
#include<iostream>
#include <eigen3/Eigen/Core>
typedef double (*pfescalarN)(Eigen::Vector2d, Eigen::Vector2d);
//Funcion eescalar -Dependiente de dos vectores- por defecto que ante cualquier input devuelve siempre 0
inline double EscalarN_Defecto(Eigen::Vector2d position,Eigen::Vector2d normal){
        return 0.0;
}
class Coeficiente_EscalarN{
    public:
    pfescalarN funcion;
    //Constructor por defecto
    Coeficiente_EscalarN(void);
    //Constructor dada una función 
    Coeficiente_EscalarN(pfescalarN FUNCION);
    //Inicializa el objeto dada una función
    virtual void Inicializa(pfescalarN FUNCION);
    //Devuelve el valor de la función cargada en punto e instante determinados
    virtual double Evalua(Eigen::Vector2d punto, Eigen::Vector2d normal);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    virtual bool Es_Numerico(void);
};
#endif