#ifndef CLASS_INTERPOLADOR
#define CLASS_INTERPOLADOR
#include "TablaDeValores.hpp"
#include <eigen3/Eigen/Core>
#include <math.h>
//Plantilla de Interpolador
class Interpolador{
    public:
    TablaDeValores valores;
    Eigen::MatrixXd matriz_interpolacion, matriz_dx, matriz_dy, matriz_dx2, matriz_dy2, matriz_dxdy,
                    matriz_dr, matriz_dt, matriz_dr2, matriz_dt2, matriz_drdt;
    //Inicializador de la clase
    Interpolador(void);
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
    //Destructor de la clase 
    ~Interpolador(void);
};
#endif