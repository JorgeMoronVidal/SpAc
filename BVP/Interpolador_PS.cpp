#include"Interpolador_PS.hpp"
void Interpolador_PS::Resetea(){
    matriz_interpolacion.resize(0,0); 
    matriz_dx.resize(0,0);  
    matriz_dy.resize(0,0);  
    matriz_dx2.resize(0,0);  
    matriz_dy2.resize(0,0);  
    matriz_dxdy.resize(0,0); 
    matriz_dr.resize(0,0); 
    matriz_dt.resize(0,0);  
    matriz_dr2.resize(0,0);  
    matriz_dt2.resize(0,0);  
    matriz_drdt.resize(0,0); 
}
void Interpolador_PS::Inicializa(TablaDeValores VALORES){
    valores = VALORES;
}
double Interpolador_PS::Interpola(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_x(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_y(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_r(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_t(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_x2(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_y2(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_r2(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_t2(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_xy(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}
double Interpolador_PS::Deriva_tr(Eigen::Vector2d punto){
    std::cout << "Programar en clases derivadas\n";
    return 0.0;
}