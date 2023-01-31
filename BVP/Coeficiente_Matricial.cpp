#include "Coeficiente_Matricial.hpp"
Coeficiente_Matricial::Coeficiente_Matricial(void){
    funcion = Matriz_Defecto;
}
Coeficiente_Matricial::Coeficiente_Matricial(pfmatriz FUNCION){
    Inicializa(FUNCION);
}
void Coeficiente_Matricial::Inicializa(pfmatriz FUNCION){
    funcion = FUNCION;
}
Eigen::Matrix2d Coeficiente_Matricial::Evalua(Eigen::Vector2d punto){
    return funcion(punto);
}
bool Coeficiente_Matricial::Es_Numerico(void){
    return false;
}