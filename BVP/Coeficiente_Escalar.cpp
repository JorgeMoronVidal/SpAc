#include "Coeficiente_Escalar.hpp"
Coeficiente_Escalar::Coeficiente_Escalar(void){
    funcion = Escalar_Defecto;
}
Coeficiente_Escalar::Coeficiente_Escalar(pfescalar FUNCION){
    Inicializa(FUNCION);
}
void Coeficiente_Escalar::Inicializa(pfescalar FUNCION){
    funcion = FUNCION;
}
double Coeficiente_Escalar::Evalua(Eigen::Vector2d punto){
    return funcion(punto);
}
bool Coeficiente_Escalar::Es_Numerico(void){
    return false;
}