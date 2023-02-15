#ifndef SISTEMALINEAL_CLASS
#define SISTEMALINEAL_CLASS
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "Interfaz.hpp"
class SistemaLineal{
    public:
    std::vector<std::vector<int>> lista_adyacencia;
    std::vector<double> diagonal;
    int n_filas,n_columnas,nnz;
    int *G_i, *G_j;
    double *G_ij, *B, *u;
    void Reset_Lista_Adyacencia(int number_knots);
    void Anadir_Lista_Adyacencia(Interfaz &interfaz);
    void Purgar_Lista_Adyacencia(void);
    void Inicializa_CSR(void);
    void Actualiza_GB(Nudo &nudo);
    void Actualiza_GB(Interfaz &interfaz);
    void Convierte_COO(std::vector<int> &G_i_COO, std::vector<int>&G_j_COO, std::vector<double> &G_ij_COO);
    void Escribe_GB_COO(void);
    void Escribe_GB_CSR(void);
    ~SistemaLineal();
};
#endif
