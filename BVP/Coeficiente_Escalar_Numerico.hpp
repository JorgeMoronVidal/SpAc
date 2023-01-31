#ifndef COEFICIENTE_ESCALAR_NUMERICO
#define COEFICIENTE_ESCALAR_NUMERICO
#include"Coeficiente_Escalar.hpp"
#include"Interpolador_PS.hpp"
#include"Interpolador_RBF.hpp"
class Coeficiente_Escalar_Numerico: public Coeficiente_Escalar{
    public:
    Interpolador interpolador;
    //Constructor por defecto
    Coeficiente_Escalar_Numerico(void);
    //Constructor dada una función 
    Coeficiente_Escalar_Numerico(Interpolador INTERPOLADOR);
    //Inicializa el objeto dada una función
    void Inicializa(Interpolador INTERPOLADOR);
    //Devuelve el valor de la función cargada en punto e instante determinados
    double Evalua(Eigen::Vector2d punto);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    bool Es_Numerico(void);
};
#endif