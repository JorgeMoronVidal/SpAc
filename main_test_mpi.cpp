#include<iostream>
#include <ctime>
#include"BVPs/Monegros_Poisson.hpp"
#include"Dominios/Circulo.hpp"
#include"Dominios/Rectangulo.hpp"
#include"Solvers/FeynmanKacSolver.hpp"
#include"Solvers/PseudoespectralCirculoSolver.hpp"
#include"Mesh/Malla.hpp"
#include"SpAc/GestorMPI.hpp"

int main(int argc, char *argv[]){
    GestorMPI gestor;
    gestor.Inicializa(argc,argv);
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
    parametros_circulo[0] = 0.0;
    parametros_circulo[1] = 0.0;
    parametros_circulo[2] = 0.0;
    std::vector<double> parametros_rectangulo;
    parametros_rectangulo.resize(4);
    parametros_rectangulo[0] = -50.0;
    parametros_rectangulo[1] = -50.0;
    parametros_rectangulo[2] = 50.0;
    parametros_rectangulo[3] = 50.0;
    bvp.frontera_dominio.Inicializa(parametros_rectangulo,Rectangulo,Stopping,Stopping);
    bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
    Nudo aux_nudo;
    Interfaz aux_interfaz;
    int N = 1000, Nr = 43, Nt = 44, Nl = 10;
    double h = 1E-03, t0 = 0.12, superposicion = 1.2;
    if(gestor.id == gestor.servidor){
        //std::cout << __FILE__ << " " << __LINE__ << " Soy el servidor\n";
        Eigen::Vector2d aux;
        Malla mesh;
        std::map<std::string,double> c2;
        c2["esquina"] = 1.5; //Cond Psi = 1.02E+10
        c2["lado"] = 1.1; //Cond Psi = 1.11E+10
        mesh.Construye_Cuadrado(parametros_rectangulo,Nl,superposicion,t0,Nt,bvp,c2);
        std::vector<Interfaz>::iterator it_interfaz = mesh.Interfaces.begin();
        int aux_indice_local = 0;
        do{
            gestor.Recibe_peticion_trabajo();
            aux_nudo.Resetea();
            aux_interfaz.Resetea();
            switch (gestor.trabajo[0])
            {
            case nuevo_trabajo:
                if((*it_interfaz).es_perimeter){
                    //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    gestor.trabajo[0] = nudo_FKAC;
                    gestor.Anuncia_trabajo();
                    //std::cout << "Destino " << gestor.status.MPI_SOURCE << std::endl;
                    gestor.Envia_interfaz((*it_interfaz),gestor.status.MPI_SOURCE);
                    //std::cout << __FILE__ << " " << __LINE__ << " " << aux_indice_local << " " << (*it_interfaz).nudos_interior.size() << std::endl;
                    gestor.Envia_nudo((*it_interfaz).nudos_interior[aux_indice_local],
                    gestor.status.MPI_SOURCE, int_nudo_FKAC, double_nudo_FKAC);
                    //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    aux_indice_local ++;
                    if(aux_indice_local == (*it_interfaz).nudos_interior.size()){
                        it_interfaz ++;
                        aux_indice_local = 0;
                    } 
                }else{
                    //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    gestor.trabajo[0] = interfaz_pseudoespectral;
                    gestor.Anuncia_trabajo();
                    gestor.Envia_interfaz((*it_interfaz),gestor.status.MPI_SOURCE);
                    it_interfaz ++;
                }
                break;
            case nudo_resuelto:
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                gestor.Recibe_nudo(aux_nudo,int_nudo_FKAC, double_nudo_FKAC);
                mesh.sistema.Actualiza_GB(aux_nudo);
                break;
            case interfaz_resuelta:
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                gestor.Recibe_interfaz(aux_interfaz);
                mesh.sistema.Actualiza_GB(aux_interfaz);
                break;
            default:
                std::cout << __FILE__ << " " << __LINE__  << " " << gestor.trabajo[0]<< " ERROR\n";
                break;
            }
        }while(it_interfaz != mesh.Interfaces.end());
        int nodes_finished = 0;
        do{
            gestor.Recibe_peticion_trabajo();
            aux_nudo.Resetea();
            aux_interfaz.Resetea();
            switch (gestor.trabajo[0]){
            case nuevo_trabajo:
                    //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    gestor.trabajo[0] = termina_construccion_G_B;
                    gestor.Anuncia_trabajo();
                    nodes_finished ++;
                break;
            case nudo_resuelto:
                gestor.Recibe_nudo(aux_nudo,int_nudo_FKAC, double_nudo_FKAC);
                mesh.sistema.Actualiza_GB(aux_nudo);
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                break;
            case interfaz_resuelta:
                gestor.Recibe_interfaz(aux_interfaz);
                mesh.sistema.Actualiza_GB(aux_interfaz);
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                break;
            default:
                //std::cout << __FILE__ << " " << __LINE__ << gestor.trabajo[0] <<  " ERROR\n";
                exit(1);
                break;
            }
        }while(nodes_finished != gestor.servidor);
        mesh.sistema.Escribe_GB_COO();
        mesh.Escribe_Posiciones(bvp);
    }else{
        //std::cout << __FILE__ << " " << __LINE__ << " Soy el trabajador " << gestor.id << std::endl;
        PseudoespectralCirculoSolver solver_ps;
        FeynmanKacSolver solver_fk(gestor.id + 1);
        do{
            gestor.trabajo[0] = nuevo_trabajo;
            gestor.Pide_trabajo();
            gestor.Recibe_anuncio_trabajo();
            aux_nudo.Resetea();
            aux_interfaz.Resetea();
            switch (gestor.trabajo[0])
            {
            case nudo_FKAC:
                gestor.Recibe_interfaz(aux_interfaz);
                gestor.Recibe_nudo(aux_nudo,int_nudo_FKAC, double_nudo_FKAC);
                parametros_circulo[0] = aux_interfaz.radio;
                parametros_circulo[1] = aux_interfaz.centro(0);
                parametros_circulo[2] = aux_interfaz.centro(1);
                bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
                solver_fk.Solve_OMP_Analytic(aux_nudo,N, h,1.0,bvp, aux_interfaz);
                gestor.trabajo[0] = nudo_resuelto;
                gestor.Pide_trabajo();
                gestor.Envia_nudo(aux_nudo, gestor.servidor, int_nudo_FKAC, double_nudo_FKAC);
                break;
            case interfaz_pseudoespectral:
                //std::cout << "Resuelve interfaz";
                gestor.Recibe_interfaz(aux_interfaz);
                solver_ps.Inicializa(aux_interfaz.centro,aux_interfaz.radio,
                aux_interfaz.nudos_circunferencia.size(),Nr,
                aux_interfaz.nudos_circunferencia[0].posicion_angular
                ,bvp);
                solver_ps.Resuelve(bvp,aux_interfaz.nudos_interior,
                aux_interfaz.nudos_circunferencia);
                gestor.trabajo[0] = interfaz_resuelta;
                gestor.Pide_trabajo();
                gestor.Envia_interfaz(aux_interfaz,gestor.servidor);
                break;
            case termina_construccion_G_B:
                break;
            default:
                std::cout << __FILE__ << " " << __LINE__ << " ERROR\n";
                exit(1);
                break;
            }
    
        }while(gestor.trabajo[0] != termina_construccion_G_B);
    }
    gestor.Finaliza();
}