#ifndef MALLA_CLASS
#define MALLA_CLASS
#include "Interfaz.hpp"
#include "Subdominio.hpp"
#include "SistemaLineal.hpp"
#include <fstream>
class Malla{
    public:
    std::vector<Interfaz> Interfaces;
    std::vector<Subdominio> Subdominios;
    std::vector<double> parametros_dominio;
    double radio, angulo_inicial, superposicion, distancia;
    int nudos_por_circulo, numero_subdominios_lado, numero_total_nudos;
    std::map<std::string, double> c2;
    std::map<std::string,int> vecinos;
    SistemaLineal sistema;
    Malla(void);
    void Reinicia(void);
    void Construye_Cuadrado(std::vector<double> PARAMETROS_DOMINIO, int NUMERO_SUBDOMINIOS_LADO, double SUPERPOSICION, 
            double ANGULO_INICIAL, int NUDOS_POR_CIRCULO, BVP problema, std::map<std::string, double> C2);
    bool Es_vecino(std::vector<Interfaz>::iterator it_int, std::map<std::string, int>::iterator it_vecinos);
    void Escribe_Posiciones(BVP problema);
    ~Malla();
};
#endif