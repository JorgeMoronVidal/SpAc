#ifndef COEFICIENTE_VECTORIAL_NUMERICO
#define COEFICIENTE_VECTORIAL_NUMERICO
#include"Coeficiente_Vectorial.hpp"
#include"Interpolador_PS.hpp"
#include"Interpolador_RBF.hpp"
class Coeficiente_Vectorial_Numerico: public Coeficiente_Vectorial{
    public:
    Interpolador interpolador_x, interpolador_y;
    Eigen::Vector2d aux;
    //Constructor por defecto
    Coeficiente_Vectorial_Numerico(void);
    //Constructor dada una función 
    Coeficiente_Vectorial_Numerico(Interpolador INTERPOLADOR_X, Interpolador INTERPOLADOR_Y);
    //Inicializa el objeto dada una función
    void Inicializa(Interpolador INTERPOLADOR_X, Interpolador INTERPOLADOR_Y);
    //Devuelve el valor de la función cargada en punto e instante determinados
    Eigen::Vector2d Evalua(Eigen::Vector2d punto);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    bool Es_Numerico(void);
};
#endif