#ifndef CIRCLE_DOMAIN
#define CIRCLE_DOMAIN
#include <vector>
#include <math.h>
#include <eigen3/Eigen/Core>
#include <algorithm>
inline double Circulo(std::vector<double> parametros, 
            Eigen::Vector2d posicion, 
            Eigen::Vector2d & exitpoint,
            Eigen::Vector2d & normal){
    
    //Give the distance to a sphere of radious parameters[0] centered in 
    //(parameters[1],parameters[2])
    normal(0) = posicion(0) - parametros[1];
    normal(1) = posicion(1) - parametros[2];
    double r = normal.norm();
    r = std::max(5E-4,r);
    normal = normal/r;
    exitpoint(0) = normal(0) * parametros[0] + parametros[1];
    exitpoint(1) = normal(1) * parametros[0] + parametros[2];
    return r - parametros[0];
}
inline bool Stopping_c(Eigen::VectorXd posicion){
    return true;
}
void Plot_Circle(double * parametros, Eigen::VectorXd & posicion){
    FILE *f;
    double N = 200.0;
    f = fopen("Output/trayectories/surface.txt", "w");
    fprintf(f,"X,Y\n");
    for(double phi = 0 ; phi <= 2 * M_PI+0.001 ; phi += (2.0/N) * M_PI){
        fprintf(f,"%e,%e\n",parametros[1] + parametros[0]*cos(phi), parametros[2] + parametros[0]*sin(phi));
    }
    fclose(f);
}
#endif
