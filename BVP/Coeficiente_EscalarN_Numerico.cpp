#include "Coeficiente_EscalarN_Numerico.hpp"
Coeficiente_EscalarN_Numerico::Coeficiente_EscalarN_Numerico(void){
}
Coeficiente_EscalarN_Numerico::Coeficiente_EscalarN_Numerico(Interpolador INTERPOLADOR){
    Inicializa(INTERPOLADOR);
}
void Coeficiente_EscalarN_Numerico::Inicializa(Interpolador INTERPOLADOR){
    interpolador = INTERPOLADOR;
}
double Coeficiente_EscalarN_Numerico::Evalua(Eigen::Vector2d punto, Eigen::Vector2d normal){
    return interpolador.Interpola(punto);
}
bool Coeficiente_EscalarN_Numerico::Es_Numerico(void){
    return true;
}