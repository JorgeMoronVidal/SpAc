#include"Coeficiente_Matricial_Numerico.hpp"
Coeficiente_Matricial_Numerico::Coeficiente_Matricial_Numerico(void){
    aux(0,0) = sqrt(-1); 
    aux(0,1) = sqrt(-1);
    aux(1,0) = sqrt(-1); 
    aux(1,1) = sqrt(-1);
}
Coeficiente_Matricial_Numerico::Coeficiente_Matricial_Numerico(Interpolador INTERPOLADOR_11, Interpolador INTERPOLADOR_12, 
                                                               Interpolador INTERPOLADOR_21, Interpolador INTERPOLADOR_22){
    Inicializa(INTERPOLADOR_11, INTERPOLADOR_12, INTERPOLADOR_21, INTERPOLADOR_22);
}
void Coeficiente_Matricial_Numerico::Inicializa(Interpolador INTERPOLADOR_11, Interpolador INTERPOLADOR_12, 
                                                Interpolador INTERPOLADOR_21, Interpolador INTERPOLADOR_22){
    interpolador_11 = INTERPOLADOR_11;
    interpolador_12 = INTERPOLADOR_12;
    interpolador_21 = INTERPOLADOR_21;
    interpolador_22 = INTERPOLADOR_22;
}
Eigen::Matrix2d Coeficiente_Matricial_Numerico::Evalua(Eigen::Vector2d punto){
    aux(0,0) = interpolador_11.Interpola(punto);
    aux(0,1) = interpolador_12.Interpola(punto);
    aux(1,0) = interpolador_21.Interpola(punto);
    aux(1,1) = interpolador_22.Interpola(punto);
    return aux;
}
bool Coeficiente_Matricial_Numerico::Es_Numerico(void){
    return true;
}