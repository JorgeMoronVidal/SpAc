#ifndef NUDO_PERIMETRAL_CLASS
#define NUDO_PERIMETRAL_CLASS
#include "Nudo.hpp"
class Nudo_Perimetral: public Nudo{
    public:
    Eigen::MatrixXd RBF_matrix;
    double c2;
    void Inicializa(void){
        std::cout << "Nudo Perimetral. To code \n" << std::endl;
    };
};
#endif