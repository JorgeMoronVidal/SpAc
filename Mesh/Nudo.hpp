#ifndef NUDO_CLASS
#define NUDO_CLASS
#include <iostream>
#include <eigen3/Eigen/Core>
#include <vector>
#include <map>
inline double Atan0(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux >= 0) return aux;
    return aux + 2.0*M_PI;
}
inline double Atan90(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux <= 0.5*M_PI) return 2.0*M_PI + aux;
    return aux;
    
}
inline double Atan180(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    return atan2(point(1),point(0));
}
inline double Atan270(Eigen::VectorXd point){
    //Devuelve el valor de la arcotangente de point(1)/point(0)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux <= -0.5*M_PI) return 2*M_PI +  aux;
    return aux;
}
class Nudo{
    public:
    //True si el nudo está sobre la frontera, false en cualquier otro caso
    bool frontera;
    //El indice se identifica con la fila de la matriz[vector] que contiene sus elementos de G[B]
    int indice_global;
    //Indice local dentro de su interfaz
    int indice_local;
    //Interfaz a la que pertenece
    std::vector<int> indice_interfaz;
    //Subdominios a los que pertenece
    std::vector<std::vector<int>> indice_subdominios;
    //Posiciones cartesiana, angular y dentro de la matriz.
    Eigen::Vector2d posicion_cartesiana;
    double posicion_angular;
    //Coeficiente B[i], valor de la solucion sobre el nodo, de la varianza, coeficiente de pearson entre phi y xi, 
    //error estadístico cometido y bias.
    double B,solucion, error, varianza, coeficiente_Pearson, error_estadistico, bias, c2;
    //Componentes de G (Formato COO)
    std::map<int,double> G;
    Nudo(void){
        
    }
    virtual void Inicializa(void){
        std::cout << " To code \n" << std::endl;
    };
    void Resetea(void){
        indice_global = 0;
        indice_local=  0;
        indice_interfaz.resize(0);
        indice_subdominios.resize(0);
        posicion_cartesiana[0] = 0.0; posicion_cartesiana[1] = 0.0;
        posicion_angular = 0.0;
        B = 0.0;
        solucion = 0.0;
        error = 0.0;
        varianza = 0.0;
        coeficiente_Pearson = 0.0;
        error_estadistico = 0.0;
        bias = 0.0;
        c2 = 0.0;
        G.clear();
    }
    void Empaca(std::vector<int> &Nudo_int, std::vector<double> &Nudo_double){
        int aux_index;
        Nudo_int.resize(7+2*indice_subdominios.size()
        +G.size());
        Nudo_int[0] = frontera ? 1:0;
        Nudo_int[1] = indice_global;
        Nudo_int[2] = indice_local;
        Nudo_int[3] = indice_interfaz[0];
        Nudo_int[4] = indice_interfaz[1];
        Nudo_int[5] = indice_subdominios.size();
        aux_index = 6;
        for(int i = 0; i < Nudo_int[5]; i ++){
            Nudo_int[aux_index] = indice_subdominios[i][0];
            Nudo_int[aux_index+1] = indice_subdominios[i][1];
            aux_index += 2;
        }
        Nudo_int[aux_index] = G.size();
        aux_index ++;
        for(std::map<int,double>::iterator it_G = G.begin(); it_G != G.end(); it_G++){
            Nudo_int[aux_index] = it_G->first;
            aux_index ++;
        }
        if(aux_index != (Nudo_int.size())) std::cout << __FILE__ << " " << __LINE__ << " ERROR\n";
        Nudo_double.resize(11+G.size());
        Nudo_double[0] = posicion_cartesiana[0];
        Nudo_double[1] = posicion_cartesiana[1];
        Nudo_double[2] = posicion_angular;
        Nudo_double[3] = B;
        Nudo_double[4] = solucion;
        Nudo_double[5] = error;
        Nudo_double[6] = varianza;
        Nudo_double[7] = coeficiente_Pearson;
        Nudo_double[8] = error_estadistico;
        Nudo_double[9] = bias;
        Nudo_double[10] = c2;
        aux_index = 11;
        for(std::map<int,double>::iterator it_G = G.begin(); it_G != G.end(); it_G++){
            Nudo_double[aux_index] = it_G->second;
            aux_index ++;
        }
        if(aux_index != (Nudo_double.size())) std::cout << __FILE__ << " " << __LINE__ << " ERROR\n";
    }
    void Desempaca(std::vector<int> &Nudo_int, std::vector<double> &Nudo_double){
        int int_index, double_index, size_G;
        frontera = Nudo_int[0] == 1 ? true:false;
        indice_global = Nudo_int[1];
        indice_local = Nudo_int[2];
        indice_interfaz.resize(2);
        indice_interfaz[0] = Nudo_int[3];
        indice_interfaz[0] = Nudo_int[4];
        indice_subdominios.resize(Nudo_int[5]);
        int_index = 6;
        for(int i = 0; i < Nudo_int[5]; i ++){
            indice_subdominios[i].resize(2);
            indice_subdominios[i][0] = Nudo_int[int_index];
            indice_subdominios[i][1] = Nudo_int[int_index+1];
            int_index += 2;
        }
        size_G = Nudo_int[int_index];
        int_index ++;
        posicion_cartesiana[0] = Nudo_double[0];
        posicion_cartesiana[1] = Nudo_double[1];
        posicion_angular = Nudo_double[2];
        B = Nudo_double[3];
        solucion = Nudo_double[4];
        error = Nudo_double[5];
        varianza = Nudo_double[6];
        coeficiente_Pearson = Nudo_double[7];
        error_estadistico = Nudo_double[8];
        bias = Nudo_double[9];
        c2 = Nudo_double[10];
        double_index = 11;
        for(int i = 0; i < size_G; i++){
            G[Nudo_int[int_index]] = Nudo_double[double_index];
            int_index ++;
            double_index ++;
        }
        if(int_index != Nudo_int.size()) std::cout << __FILE__ << " " << __LINE__ <<  " int_index "<<
        int_index << "Nudo_int.size() " << Nudo_int.size() << std::endl;
        if(double_index != Nudo_double.size()) std::cout << __FILE__ << " " << __LINE__ <<" double_index " <<
        double_index << " Nudo_double.size() " << Nudo_double.size() << std::endl; " ERROR\n";
    }
};
#endif