#include "SPDDS.hpp"
#include "poisson_3.hpp"
#include "rectangle.hpp"
int main(int argc, char *argv[]){
    Eigen::VectorXd point;
    point.resize(2);
    for(double angle = 0.0; angle < 2.0*M_PI; angle += 0.05){
        point(0) = cos(angle);
        point(1) = sin(angle);
        printf("actual angle %f Atan angle %f\n",angle, Atan(point));
    }
}
