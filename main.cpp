//g++ -g -lm -fopenmp main.cpp BVP/*.cpp Solvers/*.cpp  -march=native  -O3
#include<iostream>
#include"BVPs/Monegros_Revamped.hpp"
#include"Dominios/Circulo.hpp"
#include"Dominios/Rectangulo.hpp"
#include"Solvers/FeynmanKacSolver.hpp"
#include"Solvers/PseudoespectralCirculoSolver.hpp"
#include"Mesh/Malla.hpp"
#include <ctime>
int main(){
    Coeficiente_Escalar u,g,c,f,f_2,g_2;
    Coeficiente_Vectorial b,F;
    Coeficiente_Matricial sigma;
    u.Inicializa(Monegros_u);
    g.Inicializa(Monegros_g);
    c.Inicializa(Monegros_c);
    f.Inicializa(Monegros_f);
    sigma.Inicializa(Monegros_sigma);
    b.Inicializa(Monegros_b);
    F.Inicializa(Monegros_F);
    BVP bvp;
    bvp.u = u;
    bvp.g = g;
    bvp.c = c;
    bvp.f = f;
    bvp.sigma = sigma;
    bvp.b = b;
    bvp.F = F;
    std::vector<double> parametros_circulo;
    parametros_circulo.resize(3);
    parametros_circulo[0] = 10.0;
    parametros_circulo[1] = 0.0;
    parametros_circulo[2] = 0.0;
    std::vector<double> parametros_rectangulo;
    parametros_rectangulo.resize(4);
    parametros_rectangulo[0] = -8.0;
    parametros_rectangulo[1] = -8.0;
    parametros_rectangulo[2] = 8.0;
    parametros_rectangulo[3] = 8.0;
    bvp.frontera_dominio.Inicializa(parametros_rectangulo,Rectangulo,Stopping,Stopping);
    bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
    BVP bvp_2;
    bvp_2 = bvp;
    bvp_2.g = g_2;
    Eigen::Vector2d aux;
    Malla mesh;
    std::map<std::string,double> c2;
    mesh.Construye_Cuadrado(parametros_rectangulo,5,1.2,0.1428,44,bvp,c2);
}
   /*FeynmanKacSolver solver;
    Eigen::Vector2d X0;
    X0(0) = 1.0; X0(1) = 1.0;
    unsigned t0,t1;
    for(double h = 0.8; h > 0.01; h = h*0.5){
        t0 = clock();
        solver.Solve_OMP_Analytic(X0,10000,h,0.1,bvp);
        t1 = clock();
        printf("%f %f %f %f %f\n", bvp.u.Evalua(X0),solver.phi_VR,
        fabs(bvp.u.Evalua(X0)-solver.phi_VR), sqrt(solver.phiphi_VR- pow(solver.phi_VR,2))/sqrt(10000.0),double(t1-t0)/CLOCKS_PER_SEC );
    }*/
    /*PseudoespectralCirculoSolver PS_Solver;
    Eigen::Vector2d X0 = {1.0,2.0};
    unsigned t0, t1, t2,t3;
    t0 = clock();
    PS_Solver.Inicializa(X0,10.0,44,23,0.12831,bvp);
    t1 = clock();
    PS_Solver.Resuelve(bvp);
    t2 = clock();
    PS_Solver.Resuelve(bvp);
    t3 = clock();
    Eigen::Vector2d aux;
    for(int i = 0; i < PS_Solver.interpolador.valores.z.rows(); i++){
        for(int j = 0; j < PS_Solver.interpolador.valores.z.cols(); j++){
            aux(0) = PS_Solver.interpolador.valores.x(i,j);
            aux(1) = PS_Solver.interpolador.valores.y(i,j);
            std::cout <<  PS_Solver.interpolador.valores.z(i,j) - bvp.u.Evalua(aux)<< " ";;
        }
        std::cout << std::endl;
    }
    std::cout<< "Tiempo Inicializar " << double(t1-t0)/CLOCKS_PER_SEC << std::endl;
    std::cout<< "Tiempo Resolver " << double(t2-t1)/CLOCKS_PER_SEC << std::endl;
    std::cout<< "Tiempo Resolver " << double(t3-t2)/CLOCKS_PER_SEC << std::endl;*/