#include "Coeficiente_Escalar_Numerico.hpp"
Coeficiente_Escalar_Numerico::Coeficiente_Escalar_Numerico(void){
}
Coeficiente_Escalar_Numerico::Coeficiente_Escalar_Numerico(Interpolador INTERPOLADOR){
    Inicializa(INTERPOLADOR);
}
void Coeficiente_Escalar_Numerico::Inicializa(Interpolador INTERPOLADOR){
    interpolador = INTERPOLADOR;
}
double Coeficiente_Escalar_Numerico::Evalua(Eigen::Vector2d punto){
    return interpolador.Interpola(punto);
}
bool Coeficiente_Escalar_Numerico::Es_Numerico(void){
    return true;
}