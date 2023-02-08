#include "Malla.hpp"
Malla::Malla(void){
    Reinicia();
};

void Malla::Reinicia(void){
    Interfaces.resize(0);
    Subdominios.resize(0);
    parametros_dominio.resize(0);
    radio = 0; superposicion = 1; distancia = 0;
}
void Malla::Construye_Cuadrado(std::vector<double> PARAMETROS_DOMINIO,int NUMERO_SUBDOMINIOS_LADO, double SUPERPOSICION, 
            double ANGULO_INICIAL, int NUDOS_POR_CIRCULO, BVP problema, std::map<std::string, double> C2){
            Interfaz interfaz_centinela;
            int nudo_centinela = 0;
            parametros_dominio.resize(4);
            parametros_dominio[0] = PARAMETROS_DOMINIO[0];
            parametros_dominio[1] = PARAMETROS_DOMINIO[1];
            parametros_dominio[2] = PARAMETROS_DOMINIO[2];
            parametros_dominio[3] = PARAMETROS_DOMINIO[3];
            superposicion = SUPERPOSICION;
            angulo_inicial = ANGULO_INICIAL;
            nudos_por_circulo = NUDOS_POR_CIRCULO;
            numero_subdominios_lado = NUMERO_SUBDOMINIOS_LADO;
            for(std::map<std::string,double>::iterator it = C2.begin(); 
                it != C2.end(); it ++) c2[it->first] = it -> second;
            //disp = (parametros_dominio[2] - parametros_dominio[0])*superposition_coefficent/(2*superposition_coefficent + 2*K_EXCESS*(N_centers-1));
            //d = (parametros_dominio[2] - parametros_dominio[0] - 2*disp)/(N_centers -1);
            distancia = (parametros_dominio[2] - parametros_dominio[0])/(numero_subdominios_lado -1);
            if(superposicion > 2.0){
                printf("FATAL WARNING: Superposition coefficent has to be equal or lower to 2.0\n");
            }
            if(superposicion < 1.0){
                printf("FATAL WARNING: Superposition coefficent has to be equal or higher than 1.0\n");
            }
            radio = superposicion*distancia/sqrt(2.0);
            Eigen::Vector2d aux_vec;
            //Initialization of all the centers
            for(int x_index = 0; x_index <= numero_subdominios_lado-1; x_index++){
                for(int y_index = 0; y_index <= numero_subdominios_lado-1; y_index++){
                    aux_vec[0] = parametros_dominio[0] + x_index*distancia; aux_vec[1] = parametros_dominio[1] + y_index*distancia;
                    nudo_centinela = interfaz_centinela.Inicia_Cuadrado(radio,problema,parametros_dominio, angulo_inicial,c2,
                    aux_vec,x_index,y_index,numero_subdominios_lado-1,numero_subdominios_lado-1, nudos_por_circulo,
                    nudo_centinela);
                    interfaz_centinela.es_perimeter = (x_index*y_index == 0 || x_index == numero_subdominios_lado-1 || y_index == numero_subdominios_lado-1);
                    Interfaces.push_back(interfaz_centinela);
                }
            }
            numero_total_nudos = nudo_centinela;
            sistema.Reset_Lista_Adyacencia(numero_total_nudos);
            double d_centinel;
            vecinos["SW"] = -numero_subdominios_lado -1;
            vecinos["W"] = -numero_subdominios_lado;
            vecinos["NW"] = -numero_subdominios_lado +1;
            vecinos["S"] = -1;
            vecinos["N"] = +1;
            vecinos["SE"] = numero_subdominios_lado -1;
            vecinos["E"] = numero_subdominios_lado;
            vecinos["NE"] = numero_subdominios_lado +1;
            std::vector<int> k_centinela;
            int n = nudos_por_circulo/2;
            for(int i = -n; i < n-1; i++) k_centinela.push_back(i); 
            //The center to which each knot belongs is identified and stored in the first element of knots
            for(std::vector<Interfaz>::iterator it_inter = Interfaces.begin(); it_inter != Interfaces.end(); it_inter++){
                //The center or centers in whic each knot is integrated are identified and stored the the latter elements 
                //Inicializar theta centinela.
                for(std::map<std::string,int>::iterator it_dir = vecinos.begin();
                it_dir != vecinos.end(); it_dir ++){
                    if(Es_vecino(it_inter,it_dir)){
                        for(std::vector<Nudo>::iterator it_nudo = (*(it_inter + it_dir->second)).nudos_circunferencia.begin();
                        it_nudo != (*(it_inter + it_dir->second)).nudos_circunferencia.end();
                        it_nudo ++){
                            if(((*it_nudo).posicion_cartesiana-(*it_inter).centro).norm()<(*it_inter).radio){
                                (*it_inter).nudos_interior.push_back(*it_nudo);
                            }
                        }
                    }
                }
            sistema.Anadir_Lista_Adyacencia(*it_inter);
            }
        sistema.Purgar_Lista_Adyacencia();
        sistema.Inicializa_CSR();
}
bool Malla::Es_vecino(std::vector<Interfaz>::iterator it_inter, std::map<std::string, int>::iterator it_vecinos){
     if ((it_vecinos->first).find('N')<(it_vecinos->first).length()){
        if ((*it_inter).posicion_centro[1] == numero_subdominios_lado - 1) return false;
     }
     if ((it_vecinos->first).find('S')<(it_vecinos->first).length()){
        if ((*it_inter).posicion_centro[1] == 0) return false;
     }
     if ((it_vecinos->first).find('E')<(it_vecinos->first).length()){
        if ((*it_inter).posicion_centro[0] == numero_subdominios_lado - 1) return false;
     }
     if ((it_vecinos->first).find('W')<(it_vecinos->first).length()){
        if ((*it_inter).posicion_centro[0] == 0) return false;
     }
     return true;
}
void Malla::Escribe_Posiciones(BVP problema){
    std::ofstream file_position("knot_position.txt");
    file_position.precision(8);
    file_position.setf(std::ios::fixed, std::ios::scientific);
    for(std::vector<Interfaz>::iterator it_inter = Interfaces.begin();
    it_inter !=  Interfaces.end(); it_inter ++){
        for(std::vector<Nudo>::iterator it_nudo = (*it_inter).nudos_circunferencia.begin();
        it_nudo != (*it_inter).nudos_circunferencia.end(); it_nudo ++){
            file_position <<(*it_nudo).indice_global <<" "<< (*it_nudo).posicion_cartesiana(0) << " " << 
            (*it_nudo).posicion_cartesiana(1) <<" " << problema.u.Evalua((*it_nudo).posicion_cartesiana) << std::endl;
        }
    }
    file_position.close();
}
Malla::~Malla(void){
    Interfaces.resize(0);
    Subdominios.resize(0);
    parametros_dominio.resize(0);
    radio = 0; superposicion = 0; distancia = 0;
}