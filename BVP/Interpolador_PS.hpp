#ifndef CLASS_INTERPOLADOR_PS
#define CLASS_INTERPOLADOR_PS
#include "Interpolador.hpp"

//Plantilla de Interpolador Pseudoespectral
class Interpolador_PS : public Interpolador{
    public:
    //Variables relacionadas con la interpolacion Baricentrica y de Fourier
    int m,n,N;
    double normN;
    Eigen::VectorXd c,diametro;
    //Resetea el interpolador
    void Resetea(void);
    //Prepara el interpolador
    void Inicializa(TablaDeValores VALORES);
    //Devuelve el valor interpolado en un punto
    double Interpola(Eigen::Vector2d punto);
    //Devuelve la derivada con respecto a x en un punto
    double Deriva_x(Eigen::Vector2d punto);
    //Devuelve la derivada con respecto a y en un punto 
    double Deriva_y(Eigen::Vector2d punto);
    //Devuelve la derivada con respecto al radio en un punto
    double Deriva_r(Eigen::Vector2d punto);
    //devuelve la derivada con respecto a la coordenada angular en un punto
    double Deriva_t(Eigen::Vector2d punto);
    //Devuelve la derivada segunda con respecto a x en un punto
    double Deriva_x2(Eigen::Vector2d punto);
    //Devuelve la derivada segunda con respecto a y en un punto 
    double Deriva_y2(Eigen::Vector2d punto);
    //Devuelve al derivada segunda con respecto al radio en un punto
    double Deriva_r2(Eigen::Vector2d punto);
    //Devuelve la derivada segunda con respecto a la coordenada angular en un punto
    double Deriva_t2(Eigen::Vector2d punto);
    //Devuelve la derivada segunda con respecto a x e y en un punto 
    double Deriva_xy(Eigen::Vector2d punto);
    //Devuelve la derivada segunda con respecto a la coordenada angular y al radio en un punto
    double Deriva_tr(Eigen::Vector2d punto);
};
#endif