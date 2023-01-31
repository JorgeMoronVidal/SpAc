#ifndef COEFICIENTE_ESCALARN_NUMERICO
#define COEFICIENTE_ESCALARN_NUMERICO
#include"Coeficiente_EscalarN.hpp"
#include"Interpolador_PS.hpp"
#include"Interpolador_RBF.hpp"
class Coeficiente_EscalarN_Numerico: public Coeficiente_EscalarN{
    public:
    Interpolador interpolador;
    //Constructor por defecto
    Coeficiente_EscalarN_Numerico(void);
    //Constructor dada una función 
    Coeficiente_EscalarN_Numerico(Interpolador INTERPOLADOR);
    //Inicializa el objeto dada una función
    void Inicializa(Interpolador INTERPOLADOR);
    //Devuelve el valor de la función cargada en punto e instante determinados
    double Evalua(Eigen::Vector2d punto, Eigen::Vector2d normal);
    //Devuelve true si es un Coeficiente Escalar numerico, false si no lo es
    bool Es_Numerico(void);
};
#endif