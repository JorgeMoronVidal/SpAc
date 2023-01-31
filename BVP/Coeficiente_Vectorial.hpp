#ifndef COEFICIENTE_VECTORIAL
#define COEFICIENTE_VECTORIAL
#include<iostream>
#include <eigen3/Eigen/Core>
//Funcion vectorial por defecto que independientemente del imput siempre devuelve un vector de 0
typedef Eigen::Vector2d(*pfvector)(Eigen::Vector2d);
inline Eigen::Vector2d Vector_Defecto(Eigen::Vector2d position){
        return position*0.0;
}
class Coeficiente_Vectorial{
    public:
    pfvector funcion;
    Coeficiente_Vectorial(void);
    Coeficiente_Vectorial(pfvector FUNCION);
    virtual void Inicializa(pfvector FUNCION);
    virtual Eigen::Vector2d Evalua(Eigen::Vector2d punto);
    virtual bool Es_Numerico(void);
};

#endif