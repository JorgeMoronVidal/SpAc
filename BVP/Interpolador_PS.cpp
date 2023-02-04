#include"Interpolador_PS.hpp"
void Interpolador_PS::Resetea(){
    matriz_interpolacion.resize(0,0);
    m = 0;
    c.resize(0);
    diametro.resize(0);
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
    m = valores.r.cols()*2;
    N = valores.r.rows();
    n = N/2;
    normN = 1.0/N;
    c.resize(m);
    diametro.resize(m);
    diametro.head(m/2) = -1.0*Eigen::VectorXd(valores.r.row(1));
    diametro.tail(m/2) = Eigen::VectorXd(valores.r.row(1)).reverse();
    diametro = diametro/valores.radio;
    double i_ts = -1.0 +1.0/m;
    for(int i = 0; i < m; i++){
        c(i) = pow(-1,i)*sin(M_PI*(i_ts + (2.0/m)*i +1)/2);
    }
}
double Interpolador_PS::Interpola(Eigen::Vector2d punto){
    /*Diameter of the subdomain where punto lies
      Auxiliary vector for the fourier Interpolation*/
    Eigen::VectorXd u_diametro = Eigen::VectorXd::Zero(m), 
                    aux_Fourier = Eigen::VectorXd::Zero(m);
    double r_punto = (punto-valores.centro).norm()/valores.radio, theta_punto;
    theta_punto =  atan2(punto(1)-valores.centro(1),punto(0)-valores.centro(0));
    theta_punto =(theta_punto >= 0) ? theta_punto: theta_punto + 2.0*M_PI;
    //Fourier interpolation along the r constant areas
    //Todo esto se podria acelerar con FFT
    for(int l = 0; l < N; l++){
        for(int k =  -n; k < n; k ++){
            aux_Fourier = cos(k*(theta_punto-valores.t(l,1)))*Eigen::VectorXd(valores.z.row(l));
            u_diametro.head(m/2) = u_diametro.head(m/2) + aux_Fourier;
            u_diametro.tail(m/2) = u_diametro.tail(m/2) + cos(k*M_PI)*aux_Fourier;
        }
    }
    u_diametro.tail(m/2) = u_diametro.tail(m/2).reverse().eval();
    u_diametro =normN*u_diametro.reverse().eval();
    //Barycentric interpolation
    double numerador = 0, denominador = 0, temp;
    for(int i = 0; i < m; i++){
        temp = c(i)/(r_punto-diametro(i));
        numerador = numerador + temp*u_diametro(i);
        denominador = denominador + temp;
    }
    return numerador/denominador;
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