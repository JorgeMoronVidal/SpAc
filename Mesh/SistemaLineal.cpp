#include "SistemaLineal.hpp"
void SistemaLineal::Reset_Lista_Adyacencia(int numero_nudos){
    lista_adyacencia.resize(numero_nudos);
    diagonal.resize(numero_nudos,0.0);
    n_filas = numero_nudos;
}

void SistemaLineal::Anadir_Lista_Adyacencia(Interfaz &interfaz){
    for(std::vector<Nudo>::iterator it_interior = interfaz.nudos_interior.begin();
    it_interior != interfaz.nudos_interior.end(); it_interior ++){
        for(std::vector<Nudo>::iterator it_circunferencia = interfaz.nudos_circunferencia.begin();
        it_circunferencia != interfaz.nudos_circunferencia.end(); it_circunferencia ++){
            lista_adyacencia[(*it_interior).indice_global].push_back((*it_circunferencia).indice_global);
            lista_adyacencia[(*it_circunferencia).indice_global].push_back((*it_interior).indice_global);
        }
        diagonal[(*it_interior).indice_global] += 1.0;
    }
}
void SistemaLineal::Purgar_Lista_Adyacencia(void){
    nnz = 0;
    int i = 0;
    std::set<int> aux_set;
    for(std::vector<std::vector<int>>::iterator it_lista = lista_adyacencia.begin();
    it_lista != lista_adyacencia.end(); it_lista ++){
        aux_set.clear();
        for(std::vector<int>::iterator it_row = (*it_lista).begin(); 
        it_row != (*it_lista).end(); it_row ++){
            aux_set.insert(*it_row);
        }
        aux_set.insert(i);
        i++;
        (*it_lista).assign(aux_set.begin(),aux_set.end());
        std::sort((*it_lista).begin(),(*it_lista).end());
        nnz += (*it_lista).size();
    }
    //std::cout << __FILE__ << " " << __LINE__ << " " << nnz  << std::endl;
}

void SistemaLineal::Inicializa_CSR(void){
    G_i = new int[n_filas + 1];
    G_j = new int[nnz];
    G_ij = new double[nnz];
    B = new double[n_filas];
    u = new double[n_filas];
    int counter_j = 0;
    for(int fila = 0; fila < n_filas; fila++){
        B[fila] = 0.0;
        u[fila] = 0.0;
        G_i[fila] = counter_j;
        for(int columna_nz = 0; columna_nz < lista_adyacencia[fila].size(); columna_nz++){
            G_j[counter_j] = lista_adyacencia[fila][columna_nz];
            if(lista_adyacencia[fila][columna_nz] == fila){
                G_ij[counter_j] = diagonal[fila];
            }else{
                G_ij[counter_j] = 0.0;
            }
            counter_j ++;
        }
    }
    G_i[n_filas] = counter_j;
    if(G_i[n_filas]!=nnz){
        std::cout << __FILE__ << " " << __LINE__<<" "<< "ERROR "
        << counter_j << " " << nnz << std::endl;
        exit(1);
    }
}

void SistemaLineal::Actualiza_GB(Nudo &nudo){
    int j;
    B[nudo.indice_global] += nudo.B;
    //std::cout <<  __FILE__ <<" " << __LINE__<<" " << nudo.indice_global << std::endl;
    if(nudo.G.size()>0){
        for(j = G_i[nudo.indice_global]; j <G_i[nudo.indice_global + 1]; j++){
            if(G_j[j] == nudo.G.begin()->first){
                break;
            }
        }
        if(j == G_i[nudo.indice_global + 1]){
            std::cout << __FILE__ << " " << __LINE__ << " " << nudo.indice_global << " ERROR\n";
        }else{
            for(std::map<int,double>::iterator it_G = nudo.G.begin();
            it_G != nudo.G.end(); it_G++){
                G_ij[j] += it_G->second;
                j ++;
            }
        }
    }
}

void SistemaLineal::Actualiza_GB(Interfaz &interfaz){
    for(std::vector<Nudo>::iterator it_nudo = interfaz.nudos_interior.begin(); 
    it_nudo!=interfaz.nudos_interior.end(); it_nudo ++){
        Actualiza_GB(*it_nudo);
    }
}
void SistemaLineal::Escribe_GB_COO(void){
    std::ofstream file_G("G.txt"), file_B("B.txt");
    file_G.precision(8);
    file_G.setf(std::ios::fixed, std::ios::scientific);
    file_B.precision(8);
    file_B.setf(std::ios::fixed, std::ios::scientific);
    for(int fila = 0; fila < n_filas; fila ++){
        file_B << B[fila] << std::endl;
        for(int j= G_i[fila]; j< G_i[fila+1]; j++){
            file_G << fila << " " << G_j[j] << " " << G_ij[j] << std::endl;
        }
    }
    file_B.close();
    file_G.close();
}
void SistemaLineal::Escribe_GB_CSR(void){
    std::ofstream file_Gi("G_i.txt"), file_Gj("G_j.txt"), file_Gij("G_ij.txt"),file_B("B.txt");
    for(int fila = 0; fila < n_filas; fila ++){
        file_Gi << G_i[fila] << std::endl;
        file_B <<B[fila] << std::endl;
        for(int j= G_i[fila]; j< G_i[fila+1]; j++){
            file_Gj << G_j[j] << std::endl;
            file_Gij << G_ij[j] << std::endl;
        }
    }
}
SistemaLineal::~SistemaLineal(){
    delete G_i,G_j,G_ij,B,u;
}