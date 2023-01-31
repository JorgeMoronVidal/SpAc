#include"TablaDeValores.hpp"
TablaDeValores::TablaDeValores(void){
    centro[0] = 0.0; centro[1] = 0.0;
    radio = 0;
    x.resize(0,0);
    y.resize(0,0);
    t.resize(0,0);
    r.resize(0,0);
    z.resize(0,0);
}
TablaDeValores::TablaDeValores(Eigen::Vector2d CENTRO, double RADIO, Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd R, Eigen::MatrixXd T, Eigen::MatrixXd Z){
    centro[0] = 0.0; centro[1] = 0.0;
    radio = 0;
    x.resize(0,0);
    y.resize(0,0);
    t.resize(0,0);
    r.resize(0,0);
    z.resize(0,0);
    Inicializa_TablaDeValores(CENTRO, RADIO, X, Y, R, T, Z);
}
//Inicializa la TablaDeValores con los elementos dados.
void TablaDeValores::Inicializa_TablaDeValores(Eigen::Vector2d CENTRO, double RADIO, Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd R, Eigen::MatrixXd T, Eigen::MatrixXd Z){   
    centro[0] = CENTRO[0];
    centro[1] = CENTRO[1];
    radio = RADIO;
    x = X;
    y = Y;
    r = R;
    t = T;
    z = Z;
}
TablaDeValores::~TablaDeValores(void){
    x.resize(0,0);
    y.resize(0,0);
    r.resize(0,0);
    t.resize(0,0);
    z.resize(0,0);
}