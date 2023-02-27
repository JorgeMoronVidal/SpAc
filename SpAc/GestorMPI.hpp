#ifndef CLASEGESTORMPI
#define CLASEGESTORMPI
#include <mpi.h>
#include "../Mesh/Interfaz.hpp"
#include "../Mesh/SistemaLineal.hpp"
enum EtiquetasMPI{
    peticion_trabajo,
    anuncio_trabajo,
    int_nudo_FKAC,
    double_nudo_FKAC,
    int_interfaz,
    double_interfaz,
    int_sistema,
    double_sistema,
    int_nudo_interfaz = 1000,
    double_nudo_interfaz = 1001
};
enum Trabajos{
    nuevo_trabajo,
    nudo_FKAC,
    interfaz_pseudoespectral,
    nudo_resuelto,
    interfaz_resuelta,
    bloque_G_LU,
    subdominio_pseudoespectral,
    termina_construccion_G_B,
    termina_GMRES
};
namespace MPIUtils{
    template <typename T>
    MPI_Datatype typeMPI(void);
    template <>
    MPI_Datatype typeMPI<double>(void){
        return MPI_DOUBLE;
    }
    template <>
    MPI_Datatype typeMPI<int>(void){
        return MPI_INT;
    }
}
class GestorMPI{
    public:
    MPI_Status status;
    int trabajo[1],numprocesos, id, servidor, tipo_multithreading;
    GestorMPI(){};
    GestorMPI(int argc, char *argv[]){
        Inicializa(argc,argv);
    };
    void Inicializa(int argc, char *argv[]){
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,& tipo_multithreading);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocesos);
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        servidor = numprocesos-1;
    };
    template <typename T>
    inline MPI_Status Envia_vector(std::vector<T> & vector, int destino, int etiqueta){
        MPI_Send(&vector[0],vector.size(),MPIUtils::typeMPI<T>(),destino,etiqueta,MPI_COMM_WORLD);
        return status;
    };
    template<typename T>
    inline MPI_Status Recibe_vector(std::vector<T> & vector, int etiqueta){
        MPI_Probe(MPI_ANY_SOURCE, etiqueta, MPI_COMM_WORLD, &status);
        int numero_elementos;
        MPI_Get_count(&status, MPIUtils::typeMPI<T>(), &numero_elementos);
        vector.resize(numero_elementos);
        MPI_Recv(&vector[0],vector.size(),MPIUtils::typeMPI<T>(),status.MPI_SOURCE,etiqueta,MPI_COMM_WORLD,& status);
        return status;
    };
    template<typename T>
    inline MPI_Status Recibe_vector(std::vector<T> & vector, int origen, int etiqueta){
        MPI_Probe(origen, etiqueta, MPI_COMM_WORLD, &status);
        int numero_elementos;
        MPI_Get_count(&status, MPIUtils::typeMPI<T>(), &numero_elementos);
        vector.resize(numero_elementos);
        MPI_Recv(&vector[0],vector.size(),MPIUtils::typeMPI<T>(),origen,etiqueta,MPI_COMM_WORLD,& status);
        return status;
    };
    template<typename T>
    void Retransmite_vector(std::vector<T> & vector, int origen){
       MPI_Bcast(&vector[0],vector.size(), MPIUtils::typeMPI<T>(),origen, MPI_COMM_WORLD);
    };
    inline MPI_Status Pide_trabajo(){
        MPI_Send(trabajo,1,MPI_INT,servidor,peticion_trabajo,MPI_COMM_WORLD);
        return status;
    };
    inline void Recibe_peticion_trabajo(){
        MPI_Recv(trabajo,1,MPI_INT,MPI_ANY_SOURCE,peticion_trabajo,MPI_COMM_WORLD,&status);
    };
    inline MPI_Status Anuncia_trabajo(){
        MPI_Send(trabajo,1,MPI_INT,status.MPI_SOURCE,anuncio_trabajo,MPI_COMM_WORLD);
        return status;
    };
    inline void Recibe_anuncio_trabajo(){
        MPI_Recv(trabajo,1,MPI_INT,MPI_ANY_SOURCE,anuncio_trabajo,MPI_COMM_WORLD,&status);
    };
    inline MPI_Status Envia_nudo(Nudo & nudo, int destino, int etiqueta_int, int etiqueta_double){
        std::vector<int> Nudo_int;
        std::vector<double> Nudo_double;
        nudo.Empaca(Nudo_int,Nudo_double);
        Envia_vector<int>(Nudo_int,destino,etiqueta_int);
        Envia_vector<double>(Nudo_double,destino,etiqueta_double);
        Nudo_int.clear();
        Nudo_double.clear();
        return status;
    };
    inline MPI_Status Recibe_nudo(Nudo & nudo, int origen, int etiqueta_int, int etiqueta_double){
        std::vector<int> Nudo_int;
        std::vector<double> Nudo_double;
        Recibe_vector<int>(Nudo_int,origen,etiqueta_int);
        Recibe_vector<double>(Nudo_double,origen,etiqueta_double);
        nudo.Desempaca(Nudo_int, Nudo_double);
        return status;
    };
    inline MPI_Status Envia_interfaz(Interfaz & interfaz, int destino){
        std::vector<int> Interfaz_int;
        std::vector<double> Interfaz_double;
        std::vector<Nudo> Interfaz_nudo;
        interfaz.Empaca(Interfaz_int,Interfaz_double,Interfaz_nudo);
        Envia_vector<int>(Interfaz_int,destino, int_interfaz);
        Envia_vector<double>(Interfaz_double, destino, double_interfaz);
        for(int i = 0; i < Interfaz_int[7] + Interfaz_int[8]; i++){
            Envia_nudo(Interfaz_nudo[i],destino,int_nudo_interfaz + i,double_nudo_interfaz + i);
        }
        return status;
    };
    inline MPI_Status Recibe_interfaz(Interfaz & interfaz, int origen){
        std::vector<int> Interfaz_int;
        std::vector<double> Interfaz_double;
        std::vector<Nudo> Interfaz_nudo;
        Recibe_vector<int>(Interfaz_int,origen,int_interfaz);
        Recibe_vector<double>(Interfaz_double,status.MPI_SOURCE,double_interfaz);
        Interfaz_nudo.resize(Interfaz_int[7] + Interfaz_int[8]);
        for(int i = 0; i < Interfaz_int[7] + Interfaz_int[8]; i++){
            Recibe_nudo(Interfaz_nudo[i],origen,int_nudo_interfaz + i,double_nudo_interfaz +i);
        }
        interfaz.Desempaca(Interfaz_int,Interfaz_double,Interfaz_nudo);
        return status;
    };
    MPI_Status Envia_Sistema_Lineal(SistemaLineal & sistema, int destino){
        std::vector<int> sistema_int;
        std::cout << int(3+sistema.n_filas + sistema.nnz) << std::endl;
        sistema_int.resize(int(3+sistema.n_filas + sistema.nnz));
        std::vector<double> sistema_double;
        sistema_double.resize(int(sistema.n_filas + sistema.nnz));
        sistema_int[0] = sistema.n_filas;
        sistema_int[1] = sistema.nnz;
        int i;
        for(i = 0; i < sistema.n_filas; i++){
            sistema_int[2+i] = sistema.G_i[i];
            sistema_double[i] = sistema.B[i];
        }
        sistema_int[2+i] = sistema.G_i[i];
        for(i = 0; i <sistema.nnz; i++){
            sistema_int[3+sistema.n_filas+i] = sistema.G_j[i];
            sistema_double[sistema.n_filas+i] = sistema.G_ij[i];
        }
        Envia_vector<int>(sistema_int,destino,int_sistema);
        Envia_vector<double>(sistema_double,destino,double_sistema);
        return status;
    };
    MPI_Status Recibe_Sistema_Lineal(SistemaLineal & sistema, int origen){
        std::vector<int> sistema_int;
        std::vector<double> sistema_double;
        Recibe_vector<int>(sistema_int,origen,int_sistema);
        Recibe_vector<double>(sistema_double,origen,double_sistema);
        sistema.n_filas = sistema_int[0];
        sistema.nnz = sistema_int[1];
        sistema.B = new double[sistema.n_filas];
        sistema.u = new double[sistema.n_filas];
        sistema.G_i = new int[sistema.n_filas +1];
        sistema.G_j = new int[sistema.nnz];
        sistema.G_ij = new double[sistema.nnz];
        int i;
        for(i = 0; i < sistema.n_filas; i++){
            sistema.G_i[i] = sistema_int[2+i];
            sistema.B[i] = sistema_double[i];
        }
        sistema.G_i[i] = sistema_int[2+i];
        for(i = 0; i <sistema.nnz; i++){
            sistema.G_j[i] = sistema_int[3+sistema.n_filas+i];
            sistema.G_ij[i] = sistema_double[sistema.n_filas+i];
        }
        return status;
    };
    void Retransmite_Sistema_Lineal(SistemaLineal & sistema, int origen){
        int aux_int[2];
        if(id == origen){
            aux_int[0] = sistema.n_filas;
            aux_int[1] = sistema.nnz;
            MPI_Bcast(aux_int,2, MPI_INT,origen, MPI_COMM_WORLD);
            MPI_Bcast(sistema.B,sistema.n_filas, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            //MPI_Bcast(sistema.u,sistema.n_filas, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            MPI_Bcast(sistema.G_ij,sistema.nnz, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            MPI_Bcast(sistema.G_i,sistema.n_filas + 1, MPI_INT,origen, MPI_COMM_WORLD);
            MPI_Bcast(sistema.G_j,sistema.nnz, MPI_INT,origen, MPI_COMM_WORLD);
            
        }else{
            MPI_Bcast(aux_int,2, MPI_INT,origen, MPI_COMM_WORLD);
            sistema.n_filas = aux_int[0];
            sistema.nnz = aux_int[1];
            sistema.B = new double[sistema.n_filas];
            sistema.u = new double[sistema.n_filas];
            MPI_Bcast(sistema.B,sistema.n_filas, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            //MPI_Bcast(sistema.u,sistema.n_filas, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            sistema.G_ij = new double[sistema.nnz];
            MPI_Bcast(sistema.G_ij,sistema.nnz, MPI_DOUBLE,origen, MPI_COMM_WORLD);
            sistema.G_i = new int[sistema.n_filas + 1];
            MPI_Bcast(sistema.G_i,sistema.n_filas + 1, MPI_INT,origen, MPI_COMM_WORLD);
            sistema.G_j = new int[sistema.nnz];
            MPI_Bcast(sistema.G_j,sistema.nnz, MPI_INT,origen, MPI_COMM_WORLD);
        }
        
        
    };
    void Finaliza(){
        MPI_Finalize();
    };
};
#endif