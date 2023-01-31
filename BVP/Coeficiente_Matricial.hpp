#ifndef COEFICIENTE_MATRICIAL
#define COEFICIENTE_MATRICIAL
#include<iostream>
#include <eigen3/Eigen/Core>
//Funcion vectorial por defecto que independientemente del imput siempre devuelve un vector de 0
typedef Eigen::Matrix2d(*pfmatriz)(Eigen::Vector2d);
inline Eigen::Matrix2d Matriz_Defecto(Eigen::Vector2d position){
        return 0.0*Eigen::Matrix2d::Identity();
}
class Coeficiente_Matricial{
    public:
    pfmatriz funcion;
    Coeficiente_Matricial(void);
    Coeficiente_Matricial(pfmatriz FUNCION);
    virtual void Inicializa(pfmatriz FUNCION);
    virtual Eigen::Matrix2d Evalua(Eigen::Vector2d punto);
    virtual bool Es_Numerico(void);
};

#endif