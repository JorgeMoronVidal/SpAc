#include <mpi.h>
#include "../Mesh/Interfaz.hpp"
enum EtiquetasMPI{
    peticion_trabajo,
    anuncio_trabajo,
    int_nudo_FKAC,
    double_nudo_FKAC,
    int_interfaz,
    double_interfaz,
    int_nudo_interfaz,
    double_nudo_interfaz
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
        MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,& tipo_multithreading);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocesos);
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        servidor = numprocesos-1;
        if(id == servidor){

        }
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
        return status;
    };
    inline MPI_Status Recibe_nudo(Nudo & nudo, int etiqueta_int, int etiqueta_double){
        std::vector<int> Nudo_int;
        std::vector<double> Nudo_double;
        Recibe_vector<int>(Nudo_int ,etiqueta_int);
        Recibe_vector<double>(Nudo_double, etiqueta_double);
        nudo.Desempaca(Nudo_int, Nudo_double);
    };
    inline MPI_Status Envia_interfaz(Interfaz & interfaz, int destino){
        std::vector<int> Interfaz_int;
        std::vector<double> Interfaz_double;
        std::vector<Nudo> Interfaz_nudo;
        interfaz.Empaca(Interfaz_int,Interfaz_double,Interfaz_nudo);
        Envia_vector<int>(Interfaz_int,destino, int_interfaz);
        Envia_vector<double>(Interfaz_double, destino, double_interfaz);
        for(int i = 0; i < Interfaz_int[7] + Interfaz_int[8]; i++){
            Envia_nudo(Interfaz_nudo[i],destino,int_nudo_interfaz,double_nudo_interfaz);
        }
    };
    inline MPI_Status Recibe_interfaz(Interfaz & interfaz){
        std::vector<int> Interfaz_int;
        std::vector<double> Interfaz_double;
        std::vector<Nudo> Interfaz_nudo;
        Recibe_vector<int>(Interfaz_int,int_interfaz);
        Recibe_vector<double>(Interfaz_double,double_interfaz);
        Interfaz_nudo.resize(Interfaz_int[7] + Interfaz_int[8]);
        for(int i = 0; i < Interfaz_int[7] + Interfaz_int[8]; i++){
            Recibe_nudo(Interfaz_nudo[i],int_nudo_interfaz,double_nudo_interfaz);
        }
        interfaz.Desempaca(Interfaz_int,Interfaz_double,Interfaz_nudo);
    };
    void Finaliza(){
        MPI_Finalize();
    };
};
