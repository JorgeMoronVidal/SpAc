#include<iostream>
#include <ctime>
#include"BVPs/Monegros_Poisson.hpp"
#include"Dominios/Circulo.hpp"
#include"Dominios/Rectangulo.hpp"
#include"Solvers/FeynmanKacSolver.hpp"
#include"Solvers/PseudoespectralCirculoSolver.hpp"
#include"Mesh/Malla.hpp"
#include"SpAc/GestorMPI.hpp"
#include"Solvers/SistemaLinealSolverGMRES_OMP.hpp"
void LeeArchivoConfiguracion(char file[256], int & N, int &Nr, int &Nt, int &Nl, int &maxIterGMRES, int &nparicionesMETIS, int &superposicionMETIS,
                             double & h, double & theta_0, double &superposicion, double & SW, double & NE, double & tol_GMRES);
void EscribeConfiguracion(int  N, int Nr, int Nt, int Nl, int maxIterGMRES, int nparticionesMETIS, int superposicionMETIS,
                          int number_processes, int number_OMP_THREADS, double  h, double  theta_0, double superposicion, 
                          double  SW, double  NE, double  tol_GMRES);
int main(int argc, char *argv[]){
    GestorMPI gestor;
    int numero_hilos;
    gestor.Inicializa(argc,argv);
    double init_time = MPI_Wtime();
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
    int N = 100, Nr = 43, Nt = 44, Nl = 20, maxIterGMRES = 100, nparticionesMETIS = 4, superposicionMETIS = 2;
    double h = 1E-02, t0 = 0.12, superposicion = 1.2, SW = -100.0,NE = 100.0, toleranciaGMRES = 1E-08;
    char filename[256];
    sprintf(filename,"configuration");
    LeeArchivoConfiguracion(filename,N,Nr,Nt,Nl,maxIterGMRES,nparticionesMETIS,superposicionMETIS,
    h,t0,superposicion,SW,NE,toleranciaGMRES);
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
    parametros_rectangulo[0] = SW;
    parametros_rectangulo[1] = SW;
    parametros_rectangulo[2] = NE;
    parametros_rectangulo[3] = NE;
    bvp.frontera_dominio.Inicializa(parametros_rectangulo,Rectangulo,Stopping,Stopping);
    bvp.frontera_subdominio.Inicializa(parametros_circulo,Circulo,Stopping,Stopping);
    std::vector<double> time(gestor.numprocesos);
    Nudo aux_nudo;
    Interfaz aux_interfaz;
    if(gestor.id == gestor.servidor){
        std::ofstream time_file("Output/times.csv");
        EscribeConfiguracion(N,Nr,Nt,Nl,maxIterGMRES,nparticionesMETIS,superposicionMETIS,
        gestor.numprocesos,omp_get_max_threads(), h ,t0,superposicion,SW,NE,toleranciaGMRES);
        Eigen::Vector2d aux;
        Malla mesh;
        std::map<std::string,double> c2;
        c2["esquina"] = 1.5; //Cond Psi = 1.02E+10
        c2["lado"] = 1.1; //Cond Psi = 1.11E+10
        double principio = MPI_Wtime();
        time[gestor.id] = principio;
        mesh.Construye_Cuadrado(parametros_rectangulo,Nl,superposicion,t0,Nt,bvp,c2);
        time_file <<"Construye_Cuadrado,"<<gestor.id<<","<<MPI_Wtime()-time[gestor.id]<<std::endl;
        std::vector<Interfaz>::iterator it_interfaz = mesh.Interfaces.begin();
        int aux_indice_local = 0;
        do{
            gestor.Recibe_peticion_trabajo();
            aux_nudo.Resetea();
            aux_interfaz.Resetea();
            switch (gestor.trabajo[0])
            {
            case nuevo_trabajo:
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
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
                        printf("{G_ij,B_i} for Interface [%d %d] computed s\n",
                        (*it_interfaz).posicion_centro[0],(*it_interfaz).posicion_centro[1]);
                        it_interfaz ++;
                        aux_indice_local = 0;
                    } 
                }else{
                    //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    gestor.trabajo[0] = interfaz_pseudoespectral;
                    gestor.Anuncia_trabajo();
                    gestor.Envia_interfaz((*it_interfaz),gestor.status.MPI_SOURCE);
                    printf("Interface [%d %d] sent\n",
                        (*it_interfaz).posicion_centro[0],(*it_interfaz).posicion_centro[1]);
                    it_interfaz ++;
                }
                break;
            case nudo_resuelto:
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                gestor.Recibe_nudo(aux_nudo,gestor.status.MPI_SOURCE,int_nudo_FKAC, double_nudo_FKAC);
                time_file <<"Knot_FKAC,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                mesh.sistema.Actualiza_GB(aux_nudo);
                time_file <<"Update_GB,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                break;
            case interfaz_resuelta:
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                gestor.Recibe_interfaz(aux_interfaz,gestor.status.MPI_SOURCE);
                time_file <<"Interface_int,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                mesh.sistema.Actualiza_GB(aux_interfaz);
                time_file <<"Update_GB,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
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
                    time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                    gestor.Anuncia_trabajo();
                    nodes_finished ++;
                break;
            case nudo_resuelto:
                gestor.Recibe_nudo(aux_nudo,gestor.status.MPI_SOURCE,int_nudo_FKAC, double_nudo_FKAC);
                time_file <<"Knot_FKAC,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                mesh.sistema.Actualiza_GB(aux_nudo);
                time_file <<"Update_GB,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                break;
            case interfaz_resuelta:
                gestor.Recibe_interfaz(aux_interfaz, gestor.status.MPI_SOURCE);
                time_file <<"Interfaz_int,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                mesh.sistema.Actualiza_GB(aux_interfaz);
                time_file <<"Update_GB,"<<gestor.status.MPI_SOURCE<<","<<MPI_Wtime()-time[gestor.status.MPI_SOURCE]<<std::endl;
                time[gestor.status.MPI_SOURCE] = MPI_Wtime();
                //std::cout << __FILE__ << " " << __LINE__ << std::endl;
                break;
            default:
                //std::cout << __FILE__ << " " << __LINE__ << gestor.trabajo[0] <<  " ERROR\n";
                exit(1);
                break;
            }
        }while(nodes_finished != gestor.servidor);
        time_file <<"G_and_B_computation,"<<gestor.id<<","<<MPI_Wtime()-principio<<std::endl;
        time[gestor.id] = MPI_Wtime();
        //mesh.sistema.Escribe_GB_COO();
        mesh.Escribe_Posiciones(bvp);
        time_file <<"Writing_IFiles,"<<gestor.id<<","<<MPI_Wtime()-time[gestor.id]<<std::endl;
        time[gestor.id] = MPI_Wtime();
        gestor.Retransmite_Sistema_Lineal(mesh.sistema, gestor.servidor);
        time_file <<"Sending_LS,"<<gestor.id<<","<<MPI_Wtime()-time[gestor.id]<<std::endl;
        numero_hilos = omp_get_num_threads();
        MPI_Barrier(MPI_COMM_WORLD);
        //omp_set_num_threads(1);
        SistemaLinealSolverGMRES GMRESSolver;
        time[gestor.id] = MPI_Wtime();
        GMRESSolver.Resuelve(mesh.sistema,gestor,maxIterGMRES,toleranciaGMRES,nparticionesMETIS,superposicionMETIS);
        MPI_Barrier(MPI_COMM_WORLD);
        //omp_set_num_threads(numero_hilos);
        time_file <<"Solving_LS,"<<gestor.id<<","<<MPI_Wtime()-time[gestor.id]<<std::endl;
        time_file <<"Total,"<<gestor.id<<","<<MPI_Wtime()-principio<<std::endl;
        time_file.close();
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
                gestor.Recibe_interfaz(aux_interfaz,gestor.servidor);
                //std::cout << "Nudo recibido " << MPI_Wtime() - init_time << std::endl;
                gestor.Recibe_nudo(aux_nudo, gestor.servidor, int_nudo_FKAC, double_nudo_FKAC);
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
                //std::cout << "Interfaz recibida " << MPI_Wtime() - init_time << std::endl;
                gestor.Recibe_interfaz(aux_interfaz,gestor.servidor);
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
            std::cout << "Terminada construccion G_B " << gestor.id << std::endl;
                break;
            default:
                std::cout << __FILE__ << " " << __LINE__ << " ERROR\n";
                exit(1);
                break;
            }
    
        }while(gestor.trabajo[0] != termina_construccion_G_B);
        SistemaLineal sistema;
        gestor.Retransmite_Sistema_Lineal(sistema, gestor.servidor);
        //std::cout <<"Sistema Lineal Recibido " << gestor.id << std::endl;
        numero_hilos = omp_get_num_threads();
        MPI_Barrier(MPI_COMM_WORLD);
        //omp_set_num_threads(1);
        SistemaLinealSolverGMRES GMRESSolver;
        GMRESSolver.Resuelve(sistema,gestor,maxIterGMRES,toleranciaGMRES,nparticionesMETIS,superposicionMETIS);
        MPI_Barrier(MPI_COMM_WORLD);
        //omp_set_num_threads(numero_hilos);
    }
    gestor.Finaliza();
}

void LeeArchivoConfiguracion(char file[256] , int & N, int &Nr, int &Nt, int &Nl, int &maxIterGMRES, int &nparicionesMETIS, int &superposicionMETIS,
    double & h, double & theta_0, double &superposicion, double & SW, double & NE, double & tol_GMRES){
    std::ifstream myfile(file);
    std::string line, cent;
    if (myfile)  // same as: if (myfile.good())
    {   while (getline( myfile, line )){  // same as: while (getline( myfile, line ).good())
            if (line.c_str()[0]!='%') //First of all we check if the line is a comment
            { 
                cent.clear(); //Centinel needs to be empty
                for(std::string::iterator it = line.begin(); it != line.end(); it++){
                //Loop over all the characters of the line
                    if(*it != ' ' && *it != '\t'){
                    //If it is not a space nor a tabulation then we concatenate the term to cent
                        cent += *it;
                    }
                    else{
                    //If we get a space or a tabulation then we check the string
                        if(cent == "N_trayectorias="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;
                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            N = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "N_teta="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            Nt = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "N_radio="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            Nr = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "N_lado="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            Nl = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "MaxIGMRES="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            maxIterGMRES = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "NparMETIS="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;
                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            nparicionesMETIS = atoi(cent.c_str());
                            break;
                        }
                        if(cent == "NsupMETIS="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            superposicionMETIS= atoi(cent.c_str());
                            break;
                        }
                        if(cent == "h="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            h = atof(cent.c_str());
                            break;
                        }
                        if(cent == "teta_0="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            theta_0 = atof(cent.c_str());
                            break;
                        }
                        if(cent == "supsub="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            superposicion = atof(cent.c_str());
                            break;
                        }
                        if(cent == "SW="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            SW = atof(cent.c_str());
                            break;
                        }
                        if(cent == "NE="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            NE = atof(cent.c_str());
                            break;
                        }
                        if(cent == "tol_GMRES="){
                        while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            while( it != line.end()){
                                cent += *it;
                                it++;
                            }
                            tol_GMRES = atof(cent.c_str());
                            break;
                        }
                    }
                }
            }
        }
        myfile.close();
    } else {
        std::cout << "WARNING: Cofiguration file was not found\n";
        std::cout << "Make sure it exist a configuration.txt file in the program's directory.\n";
    }
}
void EscribeConfiguracion(int  N, int Nr, int Nt, int Nl, int maxIterGMRES, int nparticionesMETIS, int superposicionMETIS,
                          int number_processes, int number_OMP_THREADS, double  h, double  theta_0, double superposicion, 
                          double  SW, double  NE, double  tol_GMRES){
printf("********************************SPAC SOLVER********************************\n");
printf("**********************************ver 0.1**********************************\n");
printf("***************************** Jorge Moron-Vidal****************************\n");
printf("*****************************jmoron@math.uc3m.es***************************\n");
printf("***************************************************************************\n\n");
printf("*******************************MPI Parameters******************************\n\n");
printf("        Number of proceses = %d      Number of threads = %d   \n",
       number_processes, number_OMP_THREADS);
printf("***********************Domain Decomposition Parameters*********************\n\n");
printf("Domain:[%.1f %.1f]^2    # of subdomains = [%dX%d]    Domain superposition %.0f %%\n",
      SW,NE,Nl,Nl,100*superposicion-100);
printf("# knots per circunference: %d   # knots per radious %d  initial \theta = %.3f\n",
      Nt,Nr,theta_0);
printf("\n");
printf("****************************Feynman-Kac Parameters************************\n\n");
printf("Time discretization: %e   # trayectories = %d  \n",
      h,N);
printf("\n");
printf("*******************************GMRES Parameters***************************\n\n");
printf("Max. Iterartions: %d   # METIS patitions = %d   METIS superposition = %d \n",
      maxIterGMRES,nparticionesMETIS,superposicionMETIS);
}