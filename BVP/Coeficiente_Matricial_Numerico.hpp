#ifndef COEFICIENTE_MATRICIAL_NUMERICO
#define COEFICIENTE_MATRICIAL_NUMERICO
#include"Coeficiente_Matricial.hpp"
#include"Interpolador_PS.hpp"
#include"Interpolador_RBF.hpp"
class Coeficiente_Matricial_Numerico: public Coeficiente_Matricial{
    public:
    Interpolador interpolador_11, interpolador_12, interpolador_21, interpolador_22;
    Eigen::Matrix2d aux;
    Coeficiente_Matricial_Numerico(void);
    Coeficiente_Matricial_Numerico(Interpolador INTERPOLADOR_11, Interpolador INTERPOLADOR_12, 
                                   Interpolador INTERPOLADOR_21, Interpolador INTERPOLADOR_22);
    virtual void Inicializa(Interpolador INTERPOLADOR_11, Interpolador INTERPOLADOR_12, 
                            Interpolador INTERPOLADOR_21, Interpolador INTERPOLADOR_22);
    virtual Eigen::Matrix2d Evalua(Eigen::Vector2d punto);
    virtual bool Es_Numerico(void);
};

#endif