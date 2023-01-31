#include "SPDDS.hpp"
#include "poisson_3.hpp"
#include "rectangle.hpp"
int main(int argc, char *argv[]){
    Eigen::VectorXd point, N, exitpoint;
    point.resize(2); N.resize(2); exitpoint.resize(2);
    double params[3] = {1.0, 1.0, 0.5}, dist;
    FILE *pf;
    pf = fopen("Debug/ExitPoints.csv", "w");
    fprintf(pf,"x,y,d\n");
    for(double x = 0.0; x < 2; x += 0.1){
        for(double y = -0.5; y < 1.5; y += 0.1){
            point(0) = x;
            point(1) = y;
            dist = Circle(params, point, exitpoint, N);
            fprintf(pf, "%f,%f,%f\n", exitpoint(0), exitpoint(1), dist);
        }
    }
    fclose(pf);
    system("python3 exitpoint_plot.py");
}
