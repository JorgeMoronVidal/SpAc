#ifndef SOLVER_PSEUDOESPECTRAL_CIRCULO
#define SOLVER_PSEUDOESPECTRAL_CIRCULO
#include "../BVP/BVP.hpp"
#include "../Mesh/Nudo.hpp"
#include <fstream>
#include <eigen3/Eigen/QR>
#include <ctime>
class PseudoespectralCirculoSolver{
    public:
        Eigen::MatrixXd L_pseudospectral;
        Eigen::FullPivHouseholderQR< Eigen::MatrixXd> QR;
        Eigen::VectorXd rhs_pseudospectral;
        Interpolador_PS interpolador;
        //void Compute_L_rhs_circle(std::vector<Eigen::VectorXd> positions, );
    void Inicializa(Eigen::Vector2d centro, double radio, unsigned int N_circunferencia, 
         unsigned int N_radio, double theta_0, BVP Problema);
    void Resuelve(BVP Problema);
    void Resuelve(BVP Problema, std::vector<Nudo> & nudos_interior, std::vector<Nudo> nudos_periferia);
};
#endif