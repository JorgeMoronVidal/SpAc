#include"Frontera.hpp"
Frontera::Frontera(void){
    Distancia = Distancia_Defecto;
    Absorbente = Tipo_Defecto;
    Neumann = Tipo_Defecto;
}
Frontera::Frontera(std::vector<double> PARAMETROS, pfdistancia DISTANCIA, pftipo ABSORBENTE, pftipo NEUMANN){
    Inicializa(PARAMETROS, DISTANCIA, ABSORBENTE, NEUMANN);
}
void Frontera::Inicializa(std::vector<double> PARAMETROS, pfdistancia DISTANCIA, pftipo ABSORBENTE, pftipo NEUMANN){
    parametros.resize(PARAMETROS.size());
    for(unsigned int i = 0; i < PARAMETROS.size(); i++) parametros[i] = PARAMETROS[i];
    Distancia = DISTANCIA;
    Absorbente = ABSORBENTE;
    Neumann = NEUMANN;
}
double Frontera::Evalua_Distancia(Eigen::Vector2d punto, Eigen::Vector2d & punto_mas_cercano, Eigen::Vector2d & normal){
    return Distancia(parametros, punto, punto_mas_cercano, normal);
}
bool Frontera::Evalua_Absorbente(Eigen::Vector2d punto){
    return Absorbente(punto);
}
bool Frontera::Evalua_Neumann(Eigen::Vector2d punto){
    return Neumann(punto);
}
bool Frontera::Es_Numerico(void){
    return false;
}