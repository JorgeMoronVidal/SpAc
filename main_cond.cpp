#include<iostream>
#include <ctime>
#include"BVPs/Monegros_Poisson.hpp"
#include"Dominios/Circulo.hpp"
#include"Dominios/Rectangulo.hpp"
#include"Solvers/FeynmanKacSolver.hpp"
#include"Solvers/PseudoespectralCirculoSolver.hpp"
#include"Mesh/Malla.hpp"
#include"Solvers/SistemaLinealSolverEIGEN.hpp"
void LeeArchivoConfiguracion(char file[256], int & N, int &Nr, int &Nt, int &Nl, int &maxIterGMRES, int &nparicionesMETIS, int &superposicionMETIS,
                             double & h, double & theta_0, double &superposicion, double & SW, double & NE, double & tol_GMRES);
int main(){
    SistemaLineal sistema;
    sistema.Lee_GB_COO();
    //sistema.Escribe_GB_COO();
    SistemaLinealSolverEIGEN solver;
    solver.Initialize(sistema);
    solver.Estima_Numero_Condicionamiento(100,30);
};
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