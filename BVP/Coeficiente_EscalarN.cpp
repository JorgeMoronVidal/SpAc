#include "Coeficiente_EscalarN.hpp"
Coeficiente_EscalarN::Coeficiente_EscalarN(void){
    funcion = EscalarN_Defecto;
}
Coeficiente_EscalarN::Coeficiente_EscalarN(pfescalarN FUNCION){
    Inicializa(FUNCION);
}
void Coeficiente_EscalarN::Inicializa(pfescalarN FUNCION){
    funcion = FUNCION;
}
double Coeficiente_EscalarN::Evalua(Eigen::Vector2d punto, Eigen::Vector2d normal){
    return funcion(punto,normal);
}
bool Coeficiente_EscalarN::Es_Numerico(void){
    return false;
}