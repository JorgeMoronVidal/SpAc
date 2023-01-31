#ifndef CLASS_TABLADEVALORES
#define CLASS_TABLADEVALORES
#include <iostream>
#include <eigen3/Eigen/Core>
class TablaDeValores{
    public:
    /*Tablas de valores de los puntos a evaluar.
    Coordenadas dadas en cartesianas x e y y polares t y r.
    Valores de la función escalar almacenados en z*/
    Eigen::MatrixXd x,y,t,r,z;
    Eigen::Vector2d centro;
    double radio;
    //Constructores de la clase
    TablaDeValores(void);
    TablaDeValores(Eigen::Vector2d CENTRO, double RADIO, Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd R, Eigen::MatrixXd T, Eigen::MatrixXd Z);
    //Inicializa la TablaDeValores con los elementos dados.
    void Inicializa_TablaDeValores(Eigen::Vector2d CENTRO, double RADIO, Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd R, Eigen::MatrixXd T, Eigen::MatrixXd Z);
    ~TablaDeValores(void);
};
#endif