#include <math.h>
#include <eigen3/Eigen/Core>
#include <algorithm>
inline double Circle(double* params, 
            Eigen::VectorXd & position, 
            Eigen::VectorXd & exitpoint,
            Eigen::VectorXd & normal){
    
    //Give the distance to a sphere of radious parameters[0] centered in 
    //(parameters[1],parameters[2])
    normal(0) = position(0) - params[1];
    normal(1) = position(1) - params[2];
    double r = normal.norm();
    r = std::max(5E-4,r);
    normal = normal/r;
    exitpoint(0) = normal(0) * params[0] + params[1];
    exitpoint(1) = normal(1) * params[0] + params[2];
    return r - params[0];
}
inline bool Stopping_c(Eigen::VectorXd position){
    return true;
}
void Plot_Circle(double * params, Eigen::VectorXd & position){
    FILE *f;
    double N = 200.0;
    f = fopen("Output/trayectories/surface.txt", "w");
    fprintf(f,"X,Y\n");
    for(double phi = 0 ; phi <= 2 * M_PI+0.001 ; phi += (2.0/N) * M_PI){
        fprintf(f,"%e,%e\n",params[1] + params[0]*cos(phi), params[2] + params[0]*sin(phi));
    }
    fclose(f);
}
