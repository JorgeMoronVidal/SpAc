#ifndef NUDO_INTERIOR_CLASS
#define NUDO_INTERIOR_CLASS
#include "Nudo.hpp"
class Nudo_Interior: public Nudo{
    public:
    //Frecuencias de la transformada discreta de Fourier
    std::vector<int> frecuencias;
    void Inicializa(void){
        std::cout << "Nudo Interior. To code \n" << std::endl;
    };
};
#endif