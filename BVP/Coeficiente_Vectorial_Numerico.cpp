#include "Coeficiente_Vectorial_Numerico.hpp"
Coeficiente_Vectorial_Numerico::Coeficiente_Vectorial_Numerico(void){
    aux(0) = sqrt(-1); aux(1) = sqrt(-1);
}
Coeficiente_Vectorial_Numerico::Coeficiente_Vectorial_Numerico(Interpolador INTERPOLADOR_X, Interpolador INTERPOLADOR_Y){
    Inicializa(INTERPOLADOR_X, INTERPOLADOR_Y);
}
void Coeficiente_Vectorial_Numerico::Inicializa(Interpolador INTERPOLADOR_X, Interpolador INTERPOLADOR_Y){
    interpolador_x = INTERPOLADOR_X;
    interpolador_y = INTERPOLADOR_Y;
}
Eigen::Vector2d Coeficiente_Vectorial_Numerico::Evalua(Eigen::Vector2d punto){
    aux(0) = interpolador_x.Interpola(punto);
    aux(1) = interpolador_y.Interpola(punto);
    return aux;
}
bool Coeficiente_Vectorial_Numerico::Es_Numerico(void){
    return true;
}