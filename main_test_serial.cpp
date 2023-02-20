//g++ -g -lm -fopenmp main.cpp BVP/*.cpp Solvers/*.cpp  -march=native  -O3
#include<iostream>
#include"BVPs/Franke_Poisson.hpp"
#include"Dominios/Circulo.hpp"
#include"Dominios/Rectangulo.hpp"
#include"Solvers/FeynmanKacSolver.hpp"
#include"Solvers/PseudoespectralCirculoSolver.hpp"
#include"Mesh/Malla.hpp"
#include"Solvers/SistemaLinealSolverMUMPS.hpp"
#include <ctime>
int main(){
    Coeficiente_Escalar u,g,c,f,f_2,g_2;
    Coeficiente_Vectorial b,F;
    Coeficiente_Matricial sigma;
    u.Inicializa(Franke_u);
    g.Inicializa(Franke_g);
    c.Inicializa(Franke_c);
    f.Inicializa(Franke_f);
    sigma.Inicializa(Franke_sigma);
    b.Inicializa(Franke_b);
    F.Inicializa(Franke_F);
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
    parametros_circulo[0] = 0.0;
    parametros_circulo[1] = 0.0;
    parametros_circulo[2] = 0.0;
    std::vector<double> parametros_rectangulo;
    parametros_rectangulo.resize(4);
    parametros_rectangulo[0] = 0.0;
    parametros_rectangulo[1] = 0.0;
    parametros_rectangulo[2] = 1.25;
    parametros_rectangulo[3] = 1.25;
    bvp.frontera_dominio.Inicializa(parametros_rectangulo,Rectangulo,Stopping,Stopping);
    bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
    Eigen::Vector2d aux;
    Malla mesh;
    std::map<std::string,double> c2;
    c2["esquina"] = 1.5; //Cond Psi = 1.02E+10
    c2["lado"] = 1.1; //Cond Psi = 1.11E+10
    mesh.Construye_Cuadrado(parametros_rectangulo,4,1.0,0.1428,44,bvp,c2);
    PseudoespectralCirculoSolver solver_ps;
    FeynmanKacSolver solver_fk;
    unsigned t0,t1,t2;
    for(std::vector<Interfaz>::iterator it_inter = mesh.Interfaces.begin();
    it_inter !=  mesh.Interfaces.end(); it_inter ++){
        if((*it_inter).es_perimeter){
            parametros_circulo[0] = (*it_inter).radio;
            parametros_circulo[1] = (*it_inter).centro(0);
            parametros_circulo[2] = (*it_inter).centro(1);
            bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
            for(std::vector<Nudo>::iterator it_nudo = (*it_inter).nudos_interior.begin();
            it_nudo != (*it_inter).nudos_interior.end(); it_nudo ++){
                t0 = clock();
                solver_fk.Solve_OMP_Analytic((*it_nudo),1000, 1E-05,1.0,bvp, (*it_inter));
                t1 = clock();
                std::cout<< "Tiempo resolviendo nudo " << (*it_nudo).indice_global << " " << 
                double(t1-t0)/CLOCKS_PER_SEC << std::endl;
                mesh.sistema.Actualiza_GB(*(it_nudo));
            }
        } else {
            t0 = clock();
            solver_ps.Inicializa((*it_inter).centro,(*it_inter).radio,
            (*it_inter).nudos_circunferencia.size(),43,
            (*it_inter).nudos_circunferencia[0].posicion_angular
            ,bvp);
            solver_ps.Resuelve(bvp,(*it_inter).nudos_interior,
            (*it_inter).nudos_circunferencia);
            t1 = clock();
            std::cout<< "Tiempo resolviendo interfaz [" << (*it_inter).posicion_centro[0]<<
            ", " <<(*it_inter).posicion_centro[1] << "] " << double(t1-t0)/CLOCKS_PER_SEC << std::endl;
            mesh.sistema.Actualiza_GB(*(it_inter));
        }
    }
    mesh.sistema.Escribe_GB_COO();
    mesh.Escribe_Posiciones(bvp);
    /*SistemaLinealSolver LUSolver;
    LUSolver.MUMPS_LU(mesh.sistema,1,1,1E-10,1);
    for(int i = 0; i < mesh.sistema.n_filas; i++){
            std::cout << LUSolver.solucion[i] << std::endl;
    }*/
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