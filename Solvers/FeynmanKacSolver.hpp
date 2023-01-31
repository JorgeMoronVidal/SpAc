#ifndef FEYNMANKACSOLVER
#define FEYNMANKACSOLVER
#include "../BVP/BVP.hpp"
#include "../Mesh/Interfaz.hpp"
#include <omp.h>
#include <boost/random.hpp>
enum sumindex
{   ScoreLinear = 0,
    ScoreSublinear = 1,
    ScoreLinearNum = 2,
    ScoreSublinearNum = 3,
    ScoreLinear2 = 4,
    ScoreSublinear2 = 5,
    ScoreLinearNum2 = 6,
    ScoreSublinearNum2 = 7,
    XiLinear = 8,
    XiSublinear = 9,
    XiLinearNum = 10,
    XiSublinearNum = 11,
    XiLinear2 = 8,
    XiSublinear2 = 9,
    XiLinearNum2 = 10,
    XiSublinearNum2 = 11,
    ScoreLinearVR = 12,
    ScoreSublinearVR = 13,
    ScoreLinearVRNum = 14,
    ScoreSublinearVRNum = 15,
    ScoreLinearVR2 = 16,
    ScoreSublinearVR2 = 17,
    ScoreLinearVRNum2 = 18,
    ScoreSublinearVRNum2 = 19,
    XiScoreLinear = 20,
    XiScoreSublinear = 21,
    XiScoreLinearNum = 22,
    XiScoreSublinearNum = 23,
    tauLinear = 24,
    tauSublinear= 25,
    RNGCalls = 26,
    AvPathLen = 27
};
enum outcome {dentro, parada_dominio, parada_subdominio, reflejada};
//Solver de problemas elípticos mediante la ecuación de Feynman-Kac
class FeynmanKacSolver{
    private:
    //Procesos y cantidades vitales para el 
    std::vector<Eigen::Vector2d> X_tau;
    std::vector<double> Y_tau, Z_tau, tau, xi_tau;
    std::vector<unsigned int> RNGcallsv;
    //Generadores de números aleatorios
    std::vector<boost::mt19937> RNG;
    std::vector<boost::normal_distribution<double>> normal_dist;
    std::vector<boost::exponential_distribution<double>> exp_dist;
    //Sumatorios de cantidades importantes para el correcto desarrollo de la simulacion
    double sums[30];
    //Actualiza el incremento browniano
    void Actualiza_Incremento(Eigen::Vector2d & increment, boost::mt19937 & rng, 
    boost::normal_distribution<double> & normalrng, double sqrth);
    //Paso temporal de X,Y,Z considerando que está cerca de condiciones de contorno absorbentes
    void Paso_Absorbente(Eigen::Vector2d & X, Eigen::Vector2d & N, double & Y, double & Z, double & xi,
    double & t, double ji_t, double h, double sqrth, Eigen::Vector2d & increment, boost::mt19937 & rng, 
    boost::normal_distribution<double> & normalrng, BVP & bvp , unsigned int & N_rngcalls);
    //Paso temporal de X,Y,Z considerando que está cerca de condiciones de contorno reflejanetes
    void Paso_Reflejante(Eigen::Vector2d & X, Eigen::Vector2d & N, Eigen::Vector2d & Npro, double & Y, double & Z ,
    double & xi, double & t, double & ji_t, double rho, double & d_k,  double h, double  sqrth, Eigen::Vector2d & increment,
    boost::mt19937 & rng, boost::normal_distribution<double> & normalrng, boost::exponential_distribution<double> & exprng
    ,BVP & bvp ,unsigned int & N_rngcalls);
    //Determina la distancia de la partícula a la frontera de su dominio de integración así como su estado con respecto al mismo
    inline outcome Dentro(double & distancia, bool & stoppingbc, bool & Neumannbc, Eigen::Vector2d & X,
    double &ji_t, double & sqrth, Eigen::Vector2d & Npro, Eigen::Vector2d & N, double Gobet_Constant, BVP & bvp);
    public:
    std::vector<int> threads;
    //Time discretization of the trajectories and initial time of the trajectory T
    double h,sqrth;
    //Number of trayectories
    unsigned int N;
    //Error objective
    double eps;
    //Important estimated quantities
    double phi, phiphi, phi_sublinear, phiphi_sublinear, phi_VR, phiphi_VR, phi_sublinearVR, phiphi_sublinearVR,
    xi, xixi, xi_sublinear, xixi_sublinear, xiphi, xiphi_sublinear,  phi_num, phiphi_num, phi_sublinear_num, 
    phiphi_sublinear_num, phi_VR_num, phiphi_VR_num, phi_sublinearVR_num, phiphi_sublinearVR_num, xi_num, xixi_num, 
    xi_sublinear_num, xixi_sublinear_num, xiphi_num, xiphi_sublinear_num,APL,RNGC;
    //Inicializador por defecto de la clase
    FeynmanKacSolver(void);
    //Inicializador de la clase cuando mezclamos MPI y OpenMP
    FeynmanKacSolver(unsigned int MPIrank);
    //Simulacion de la ecuación de Feynman-Kac usando OpenMP
    void Simulacion_OMP(Eigen::Vector2d X0, unsigned int numero_trayectorias, double discretizacion_temporal,
                   double rho, BVP Problema);
    void Reduce_Analytic(BVP Problema, unsigned int numero_trayectorias);
    void Reduce_Analytic(BVP Problema, unsigned int numero_trayectorias,Nudo nudo, Interfaz interfaz);
    void Update(void);
    void Solve_OMP_Analytic(Eigen::Vector2d X0, unsigned int numero_trayectorias, double discretizacion_temporal,
                   double rho, BVP Problema);
};
#endif