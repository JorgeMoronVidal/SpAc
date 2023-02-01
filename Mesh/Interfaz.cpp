#include "Interfaz.hpp"
int Interfaz::Inicia_Cuadrado(double RADIO, BVP problema, std::vector<double> PARAMETROS_DOMINIO, 
    double theta_0, std::map<std::string,double> c2, Eigen::Vector2d CENTRO, int indice_posicion_x,
    int indice_posicion_y, int indice_maximo_x,  int indice_maximo_y,  int numero_nudos_total, 
    int primer_nudo){
        int numero_nudos_interfaz;
        radio = RADIO;
        centro[0] = CENTRO[0]; centro[1] = CENTRO[1];
        posicion_centro.resize(2);
        posicion_centro[0] = indice_posicion_x; posicion_centro[1] = indice_posicion_y;
        parametros_dominio.resize(4);
        for(int i = 0; i < 4; i++) parametros_dominio[i] = PARAMETROS_DOMINIO[i];
        nudos_circunferencia.clear();
        nudos_interior.clear();
        double angular_increment, first_angle, last_angle;
        Eigen::Vector2d aux_vector;
        if(posicion_centro[0] == 0){
            if(posicion_centro[1] == 0){
                aux_vector(1) = parametros_dominio[1] - centro[1];
                aux_vector(0) = sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                first_angle = Atan180(aux_vector);
                aux_vector(0) = parametros_dominio[0] - centro[0];
                aux_vector(1) = sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                last_angle = Atan180(aux_vector);
                numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                //N_knots[0] = numero_nudos;
                angular_increment = 2*(M_PI)/numero_nudos_total;
                interior = -11;
            } else {
                if(posicion_centro[1] == indice_maximo_y){
                    aux_vector(0) = parametros_dominio[0] - centro[0];
                    aux_vector(1) = -sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                    first_angle = Atan180(aux_vector);
                    aux_vector(1) = parametros_dominio[3] - centro[1];
                    aux_vector(0) = sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                    last_angle = Atan180(aux_vector);
                    numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                    //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    //N_knots[0] = numero_nudos;
                    angular_increment = 2*(M_PI)/(numero_nudos_total);
                    interior = -22;
                    } else {
                        aux_vector(0) = parametros_dominio[0] - centro[0];
                        aux_vector(1) = -sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                        first_angle = Atan180(aux_vector);
                        aux_vector(1) = sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                        last_angle = Atan180(aux_vector);
                        numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                        //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        //N_knots[0] = numero_nudos;
                        angular_increment = 2*(M_PI)/(numero_nudos_total);
                        interior = -4;
                    }
                }
        } else {
            if(posicion_centro[0] == indice_maximo_x){
                if(posicion_centro[1] == 0){
                    aux_vector(0) = parametros_dominio[2] - centro[0];
                    aux_vector(1) = sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                    first_angle = Atan0(aux_vector);
                    aux_vector(1) = parametros_dominio[1] - centro[1];
                    aux_vector(0) = -sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                    last_angle = Atan0(aux_vector);
                    numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                    //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    //N_knots[0] = numero_nudos;
                    angular_increment = 2*(M_PI)/(numero_nudos_total);
                    interior = -44;
                } else {
                    if(posicion_centro[1] == indice_maximo_y){
                        aux_vector(1) = parametros_dominio[3] - centro[1];
                        aux_vector(0) = -sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                        first_angle = Atan0(aux_vector);
                        aux_vector(0) = parametros_dominio[2] - centro[0];
                        aux_vector(1) = -sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                        last_angle = Atan0(aux_vector);
                        numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                        //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        //N_knots[0] = numero_nudos;
                        angular_increment = 2*(M_PI)/(numero_nudos_total);
                        interior = -33;
                    } else {
                        aux_vector(0) = parametros_dominio[3] - centro[0];
                        aux_vector(1) = sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                        first_angle = Atan0(aux_vector);
                        aux_vector(1) = -sqrt(radio*radio - aux_vector(0)*aux_vector(0));
                        last_angle = Atan0(aux_vector);
                        numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                        //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        //N_knots[0] = numero_nudos;
                        angular_increment = 2*(M_PI)/(numero_nudos_total);
                        interior = -2;
                    }
                }
            } else {
                if(posicion_centro[1] == 0){
                    aux_vector(1) = parametros_dominio[1] - centro[1];
                    aux_vector(0) = +sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                    first_angle = Atan270(aux_vector);
                    aux_vector(0) = -sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                    last_angle = Atan270(aux_vector);
                    numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                    //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    //N_knots[0] = numero_nudos;
                    angular_increment = 2*(M_PI)/(numero_nudos_total);
                    interior = -1;
                } else {
                    if(posicion_centro[1] == indice_maximo_y){
                        aux_vector(1) = parametros_dominio[3] - centro[1];
                        aux_vector(0) = -sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                        first_angle = Atan90(aux_vector);
                        aux_vector(0) = +sqrt(radio*radio - aux_vector(1)*aux_vector(1));
                        last_angle = Atan90(aux_vector);
                        numero_nudos_interfaz = (int)(numero_nudos_total*fabs(last_angle-first_angle)/(2*M_PI));
                        //angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        //N_knots[0] = numero_nudos;
                        angular_increment = 2*(M_PI)/(numero_nudos_total);
                        interior = -3;
                        //printf("first angle = %f last angle = %f angular increment = %f\n",
                        //first_angle,last_angle,angular_increment);
                    } else {
                        nudos_circunferencia.resize(numero_nudos_total);
                        angular_increment = 2*(M_PI)/(numero_nudos_total);
                        //printf("Angular increment %f\n",angular_increment);
                        first_angle = theta_0;
                        interior = 1;
                    } 
                }
            }
        }
        //printf("Centro = [%f, %f] Posicion = [%d,%d] Interior = %d\n",centro[0], centro[1],posicion_centro[0],posicion_centro[1],interior);
        int aux_index = primer_nudo;
        if(interior > 0){
            for(int i = 0; i < nudos_circunferencia.size(); i++){
                nudos_circunferencia[i].indice_local = i;
                nudos_circunferencia[i].indice_global = aux_index;
                nudos_circunferencia[i].indice_interfaz.resize(2);
                nudos_circunferencia[i].indice_interfaz[0] = posicion_centro[0]; 
                nudos_circunferencia[i].indice_interfaz[1] = posicion_centro[1];
                (nudos_circunferencia[i]).posicion_angular = first_angle + angular_increment*i;
                //printf("%f \n",N_knots[0]*360/(2*M_PI));
                aux_vector[0] =radio* cos((nudos_circunferencia[i]).posicion_angular);
                aux_vector[1] =radio* sin((nudos_circunferencia[i]).posicion_angular);
                (nudos_circunferencia[i]).posicion_cartesiana =  centro + aux_vector;
                es_perimeter = false;
                aux_index ++;
            }
        }else{
            Nudo nudo;
            nudos_circunferencia.push_back(nudo);
            double aux_theta;
            (nudos_circunferencia[0]).posicion_angular = first_angle;
            nudos_circunferencia[0].indice_local = 0;
            nudos_circunferencia[0].indice_global = aux_index;
            nudos_circunferencia[0].indice_interfaz.resize(2);
            nudos_circunferencia[0].indice_interfaz[0] = posicion_centro[0]; 
            nudos_circunferencia[0].indice_interfaz[1] = posicion_centro[1];
             //printf("%f \n",N_knots[0]*360/(2*M_PI));
             aux_vector[0] =radio* cos((nudos_circunferencia[0]).posicion_angular);
             aux_vector[1] =radio* sin((nudos_circunferencia[0]).posicion_angular);
             (nudos_circunferencia[0]).posicion_cartesiana =  centro + aux_vector;
             aux_index ++;
            for(int i = 0; i < numero_nudos_interfaz; i++){
                    aux_theta = first_angle + theta_0 + angular_increment*i;
                //if(aux_theta > first_angle + 0.5* angular_increment && aux_theta < last_angle - 0.5* angular_increment){
                    nudos_circunferencia.push_back(nudo);
                    nudos_circunferencia[i+1].indice_local = i+1;
                    nudos_circunferencia[i+1].indice_global = aux_index;
                    nudos_circunferencia[i+1].indice_interfaz.resize(2);
                    nudos_circunferencia[i+1].indice_interfaz[0] = posicion_centro[0]; 
                    nudos_circunferencia[i+1].indice_interfaz[1] = posicion_centro[1];
                    (nudos_circunferencia[i+1]).posicion_angular = first_angle + theta_0 + angular_increment*i;
                    //printf("%f \n",N_knots[0]*360/(2*M_PI));
                    aux_vector[0] =radio* cos((nudos_circunferencia[i+1]).posicion_angular);
                    aux_vector[1] =radio* sin((nudos_circunferencia[i+1]).posicion_angular);
                    (nudos_circunferencia[i+1]).posicion_cartesiana =  centro + aux_vector;
                    aux_index ++;
                //}
            }
            nudos_circunferencia.push_back(nudo);
            (nudos_circunferencia[nudos_circunferencia.size()-1]).posicion_angular = last_angle;
            nudos_circunferencia[nudos_circunferencia.size()-1].indice_local = nudos_circunferencia.size()-1;
            nudos_circunferencia[nudos_circunferencia.size()-1].indice_global = aux_index;
            nudos_circunferencia[nudos_circunferencia.size()-1].indice_interfaz.resize(2);
            nudos_circunferencia[nudos_circunferencia.size()-1].indice_interfaz[0] = posicion_centro[0]; 
            nudos_circunferencia[nudos_circunferencia.size()-1].indice_interfaz[1] = posicion_centro[1];
            //printf("%f \n",N_knots[0]*360/(2*M_PI));
            aux_vector[0] =radio* cos((nudos_circunferencia[nudos_circunferencia.size()-1]).posicion_angular);
            aux_vector[1] =radio* sin((nudos_circunferencia[nudos_circunferencia.size()-1]).posicion_angular);
            (nudos_circunferencia[nudos_circunferencia.size()-1]).posicion_cartesiana =  centro + aux_vector;
            es_perimeter = true;
            if(interior < -10){
                Calcula_Psi(problema,c2["esquina"]);
            }else{
                Calcula_Psi(problema,c2["lado"]);
            }
            aux_index ++;
        }
        return aux_index;
}

double Interfaz::Calcula_Psi(BVP problema, double C2){
    c2 = C2;
    Psi.resize(nudos_circunferencia.size(),nudos_circunferencia.size());
    for(int i = 0; i < nudos_circunferencia.size(); i++){
        for(int j = 0; j < nudos_circunferencia.size(); j++){
            Psi(i,j) = problema.RBF(nudos_circunferencia[i].posicion_angular,nudos_circunferencia[j].posicion_angular,c2);
        }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Psi);
    iPsi = lu.inverse();
    return iPsi.norm()*Psi.norm();
}

void Interfaz::Update_G(double Y, Eigen::Vector2d X, BVP problema, std::map<int,double> &G){
    double theta;
    int i = 0;
    switch (interior)
    {
    case -11:
        theta = Atan180(X);
        break;
    case -22:
        theta = Atan180(X);
        break;
    case -33:
        theta = Atan0(X);
        break;
    case -44:
        theta = Atan0(X);
        break;
    case -1:
        theta = Atan270(X);
        break;
    case -2:
        theta = Atan0(X);
        break;
    case -3:
        theta = Atan90(X);
        break;
    case -4:
        theta = Atan180(X);
        break;
    default:
        std::cout << __FILE__ << " "<< __LINE__ <<"ERROR" << std::endl;
        break;
    }
    for(int i = 0; i < nudos_circunferencia.size(); i ++){
        for(int j = 0; j < nudos_circunferencia.size(); j++){
            G[nudos_circunferencia[i].indice_global] += 
            Y*iPsi(i,j)*problema.RBF(theta,nudos_circunferencia[i].posicion_angular,c2);
        }
        i++;
    }
} 