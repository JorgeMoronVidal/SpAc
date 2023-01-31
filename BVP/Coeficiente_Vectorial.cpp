#include "Coeficiente_Vectorial.hpp"
Coeficiente_Vectorial::Coeficiente_Vectorial(void){
    funcion = Vector_Defecto;
}
Coeficiente_Vectorial::Coeficiente_Vectorial(pfvector FUNCION){
    Inicializa(FUNCION);
}
void Coeficiente_Vectorial::Inicializa(pfvector FUNCION){
    funcion = FUNCION;
}
Eigen::Vector2d Coeficiente_Vectorial::Evalua(Eigen::Vector2d punto){
    return funcion(punto);
}
bool Coeficiente_Vectorial::Es_Numerico(void){
    return false;
}