#ifndef CENTER
#define CENTER
#include <iostream>
#include <map>
#include <string>
#include <complex>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "circle.hpp"
#include <mpi.h>
#define SUPP_CIRCLE_PARAMS 200
#define SUPP_CIRCLE_INDEX 201
#define SUPP_CIRCLE_N_KNOTS 202
#define SUPP_CIRCLE_DOUBLES 203
#define SUPP_CIRCLE_UINTS 204
#define SUPP_CIRCLE_INTS 205
#define INT_CIRCLE_PARAMS 210
#define INT_CIRCLE_INDEX 211
#define INT_CIRCLE_N_KNOTS 212
#define INT_CIRCLE_DOUBLES 213
#define INT_CIRCLE_UINTS 214
#define INT_CIRCLE_INTS 215
using namespace std::complex_literals;
inline double Atan0(Eigen::VectorXd point){
    //Returns the value of atan2 but the discontinuity is in 0 
    double aux;
    aux = atan2(point(1),point(0));
    if(aux >= 0) return aux;
    return aux + 2.0*M_PI;
}
inline double Atan90(Eigen::VectorXd point){
    //Returns the value of atan2 but in the interval[0,2\pi)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux >= 0.5*M_PI) return - 2.0*M_PI + aux;
    return aux;
    
}
inline double Atan180(Eigen::VectorXd point){
    return atan2(point(1),point(0));
}
inline double Atan270(Eigen::VectorXd point){
    //Returns the value of atan2 but in the interval[0,2\pi)
    double aux;
    aux = atan2(point(1),point(0));
    if(aux <= -0.5*M_PI) return 2*M_PI +  aux;
    return aux;
}
class Center{
    //friend class SPDDS_solver;
    friend class GMSolver;
    private:
    /*
    -parameters[0] = Radious of the circle
    -parameters[1] = x coordinate of the center
    -parameters[2] = y coordinate of the center
    */
    double parameters[3];
    //Number of knots contained in the circle
    unsigned int N_knots[2];
    //Knots contained in the circle
    std::map<unsigned int,Eigen::VectorXd> knot;
    //Frequencies of the DFT 
    std::vector<int> k;
    /*-<0 if the circle defined by this center is exterior
      ->0 if the circle defined by this center is interior
    */
    int interior;
    //B value of the knot integrated in the circle
    double BN, BBN;
    //iPsi matrix of the RBF interpolator 
    Eigen::MatrixXd iPsi;
    double c2;
    public:
    //2DIndex
    unsigned int index[2];
    //Eigen_Vector of the center 
    Eigen::VectorXd center_vec;
    //G value contained in the circle
    std::map<unsigned int,double> GN, GGN;
    //Angular position of the knots with respect to the center
    std::map<unsigned int,double> theta;
    //Default initializer
    Center(void){
        parameters[0] = parameters[1] = parameters[2] = 0.0;
        index[0] = index[1] = 0;
        N_knots[0] = 0;
    }
    unsigned int Init_Square(double R, BVP bvp, double domain_parameters[4], double fac,
        Eigen::VectorXd center, unsigned int index_x, unsigned int index_y, unsigned int number_knots, 
        unsigned int first_knot, unsigned int max_v, unsigned int max_h){
        parameters[0] = R; parameters[1] = center[0]; parameters[2] = center[1];
        index[0] = index_x; index[1] = index_y;
        double angular_increment, first_angle, last_angle;
        knot.clear();
        theta.clear();
        k.clear();
        GN.clear();
        GGN.clear();
        N_knots[1] = number_knots;
        Eigen::VectorXd aux_vector;
        aux_vector.resize(2);
        if(index_x == 0){
            if(index_y == 0){
                aux_vector(1) = domain_parameters[1] - parameters[2];
                aux_vector(0) = sqrt(R*R - aux_vector(1)*aux_vector(1));
                first_angle = Atan180(aux_vector);
                aux_vector(0) = domain_parameters[0] - parameters[1];
                aux_vector(1) = sqrt(R*R - aux_vector(0)*aux_vector(0));
                last_angle = Atan180(aux_vector);
                N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                interior = -3;
            } else {
                if(index_y == max_v){
                    aux_vector(0) = domain_parameters[0] - parameters[1];
                    aux_vector(1) = -sqrt(R*R - aux_vector(0)*aux_vector(0));
                    first_angle = Atan180(aux_vector);
                    aux_vector(1) = domain_parameters[3] - parameters[2];
                    aux_vector(0) = sqrt(R*R - aux_vector(1)*aux_vector(1));
                    last_angle = Atan180(aux_vector);
                    N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                    angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    interior = -3;
                    } else {
                        aux_vector(0) = domain_parameters[0] - parameters[1];
                        aux_vector(1) = -sqrt(R*R - aux_vector(0)*aux_vector(0));
                        first_angle = Atan180(aux_vector);
                        aux_vector(1) = sqrt(R*R - aux_vector(0)*aux_vector(0));
                        last_angle = Atan180(aux_vector);
                        N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                        angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        interior = -3;
                    }
                }
        } else {
            if(index_x == max_h){
                if(index_y == 0){
                    aux_vector(0) = domain_parameters[2] - parameters[1];
                    aux_vector(1) = sqrt(R*R - aux_vector(0)*aux_vector(0));
                    first_angle = Atan0(aux_vector);
                    aux_vector(1) = domain_parameters[1] - parameters[2];
                    aux_vector(0) = -sqrt(R*R - aux_vector(1)*aux_vector(1));
                    last_angle = Atan0(aux_vector);
                    N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                    angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    interior = -1;
                } else {
                    if(index_y == max_v){
                        aux_vector(1) = domain_parameters[3] - parameters[2];
                        aux_vector(0) = -sqrt(R*R - aux_vector(1)*aux_vector(1));
                        first_angle = Atan0(aux_vector);
                        aux_vector(0) = domain_parameters[2] - parameters[1];
                        aux_vector(1) = -sqrt(R*R - aux_vector(0)*aux_vector(0));
                        last_angle = Atan0(aux_vector);
                        N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                        angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        interior = -1;
                    } else {
                        aux_vector(0) = domain_parameters[3] - parameters[1];
                        aux_vector(1) = sqrt(R*R - aux_vector(0)*aux_vector(0));
                        first_angle = Atan0(aux_vector);
                        aux_vector(1) = -sqrt(R*R - aux_vector(0)*aux_vector(0));
                        last_angle = Atan0(aux_vector);
                        N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                        angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        interior = -1;
                    }
                }
            } else {
                //Center it is not in the corners of the domain
                if(index_y == 0){
                    aux_vector(1) = domain_parameters[1] - parameters[2];
                    aux_vector(0) = +sqrt(R*R - aux_vector(1)*aux_vector(1));
                    first_angle = Atan270(aux_vector);
                    aux_vector(0) = -sqrt(R*R - aux_vector(1)*aux_vector(1));
                    last_angle = Atan270(aux_vector);
                    N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                    angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                    interior = -4;
                } else {
                    if(index_y == max_v){
                        aux_vector(1) = domain_parameters[3] - parameters[2];
                        aux_vector(0) = +sqrt(R*R - aux_vector(1)*aux_vector(1));
                        first_angle = Atan90(aux_vector);
                        aux_vector(0) = -sqrt(R*R - aux_vector(1)*aux_vector(1));
                        last_angle = Atan90(aux_vector);
                        N_knots[0] = (int)(number_knots*fabs(last_angle-first_angle)/(2*M_PI));
                        angular_increment = (last_angle-first_angle)/(N_knots[0]-1.0);
                        interior = -2;
                    } else {
                        N_knots[0] = number_knots;
                        //We make sure N_knots[0] is odd
                        if(N_knots[0]%2 == 0) 
                        {   N_knots[0] ++;
                            N_knots[1] ++;
                        }
                        angular_increment = 2*(M_PI)/(N_knots[0]);
                        first_angle = 0.0;
                        interior = 1;
                    } 
                }
            }
        }
        printf("index = [%u %u] first a = %f last a = %f N_knots[0] = %d\n",
            index_x, index_y, first_angle, last_angle, N_knots[0]);
        unsigned int aux_index = first_knot;
        center_vec = center;
        aux_vector.resize(2);
        for(unsigned int i = 0; i < N_knots[0]; i++){
            theta[aux_index] = first_angle + angular_increment*i;
            //printf("%f \n",N_knots[0]*360/(2*M_PI));
            GN[aux_index] = 0.0;
            GGN[aux_index] = 0.0;
            aux_vector[0] =R* cos(theta[aux_index]);
            aux_vector[1] =R* sin(theta[aux_index]);
            knot[aux_index] = center + aux_vector;
            aux_index ++;
        }
        BN = 0.0;
        BBN = 0.0;
        unsigned int n;
        int k_aux;
        if(N_knots[1] % 2 == 0){
                n = (N_knots[1])/2;
                k.resize(N_knots[1]);
                k_aux = (-1)*n;
                for(unsigned int i = 0; i < N_knots[1]; i++){
                    k[i] = k_aux;
                    k_aux ++;
                }
        } else {
                n = (N_knots[1]-1)/2;
                k.resize(N_knots[1]);
                k_aux = (-1)*n;
                for(unsigned int i = 0; i < N_knots[1]; i++){
                    k[i] = k_aux;
                    k_aux ++;
                }
            }
        if(interior < 0) {
            Eigen::MatrixXd Psi, I;
            int i = 0, j = 0;
            double cond;
            c2 = pow(fac*2*M_PI/pow(N_knots[0],2),2.0);
            Psi.resize(N_knots[0],N_knots[0]);
            I.resize(N_knots[0],N_knots[0]);
            I.setIdentity();
            do{
                i = 0;
                for(std::map<unsigned int, Eigen::VectorXd>::iterator it_i = knot.begin();
                    it_i != knot.end(); it_i ++){
                        j = 0;
                    for(std::map<unsigned int, Eigen::VectorXd>::iterator it_j = knot.begin();
                    it_j != knot.end(); it_j ++){
                        Psi(i,j) = bvp.rbf.Value(it_i->second - center_vec,it_j->second - center_vec,c2);
                        if(isnan(Psi(i,j))) printf("THERE IS A NAN IN PSI\n");
                        j++;
                    }
                    i++;
                }
                Eigen::FullPivLU<Eigen::MatrixXd> lu(Psi);
                if(lu.isInvertible()){
                    iPsi = lu.inverse();
                    cond = Psi.norm()*iPsi.norm();
                    if(cond < 1E+8) c2 += 0.05;
                    if(cond > 1E+10) c2 = c2 * 0.5;
                    std::cout << c2 << std::endl;
                }else{
                    printf("FATAL ERROR: Psi Matrix is not invertible.\n");
                    iPsi = Psi*0.0;
                    c2 = 0.5*c2;
                    cond = 1E+20;
                }
            }while((cond < 1E+8) || (cond > 1E+10));
            printf("Condition number %f error %f \n",Psi.norm()*iPsi.norm(), 
                    (Psi*iPsi-I).norm());
        }
        return aux_index;
    }
    inline bool Get_knot_position(unsigned int knot_index, Eigen::VectorXd & position){
        //printf("knot_index is %u \n", knot_index);
        if (knot.find(knot_index) == knot.end()) return false;
        position = knot[knot_index];
        return true;
    }
    inline Eigen::VectorXd Get_center_position(void){
        Eigen::VectorXd vector_return;
        vector_return.resize(2);
        vector_return[0] = parameters[1];
        vector_return[1] = parameters[2];
        return vector_return;
    }
    inline double Get_radius(void){
        return parameters[0];
    }
    inline bool Is_interior(void){
        if(interior > 0) return true;
        return false;
    }
    inline bool Is_inside(Eigen::VectorXd knot_position){
        if((knot_position - Get_center_position()).norm() - parameters[0] < -1E-08) return true;
        return false;
    }
    inline double Boundary_Distance(Eigen::VectorXd knot_position){
        return (knot_position - Get_center_position()).norm() - parameters[0];
    }
    void G_Update(Eigen::VectorXd point, double Y, BVP bvp){
        double angle;
        switch(interior){
            case 1:
                angle = Atan0(point -center_vec);
                break;
            case -1:
                angle = Atan0(point - center_vec);
                break;
            case -2:
                angle = Atan90(point - center_vec);
                break;
            case -3:
                angle = Atan180(point - center_vec);
                break;
            case -4:
                angle = Atan270(point - center_vec);
                break;
            default:
                printf("Something went wrong in G_Update\n");
        }
        if(fabs(center_vec(0) + parameters[0]*cos(angle) - point(0)) + 
        fabs(center_vec(1) + parameters[0]*sin(angle) - point(1))>1E-10) 
        printf("Something went wrong in atan %d\n angle %f point [%f,%f] recovered point [%f %f]",
        interior, angle, point[0] - center_vec(0), point[1] - center_vec(1),parameters[0]*cos(angle),
        parameters[0]*sin(angle));
        /*std::complex<double> i(0.0,1.0), aux;
        for(std::map<unsigned int,double>::iterator it = GN.begin();
            it != GN.end(); it ++){
            aux = std::complex<double>(0.0,0.0);
            for(unsigned int l = 0;l < N_knots[1]; l++){
                 aux += Y*exp(k[l]*(angle-theta[it->first])*i);
            }
            GN[it->first] += (-1.0/N_knots[1])*real(aux);
            GGN[it->first] += pow((1.0/N_knots[1])*real(aux),2.0);
        }
        */
        double aux;
        if(interior == 1){
            for(std::map<unsigned int,double>::iterator it = theta.begin();
                it != theta.end(); it ++){
                aux = 0.0;
                for(unsigned int l = 0;l < N_knots[1]; l++){
                    aux += cos(k[l]*(angle-(it->second)));
                }
                GN[it->first] += Y*(-1.0/N_knots[1])*aux;
                GGN[it->first] += pow((Y/N_knots[1])*aux,2.0);
            }
        } else {
            int j = 0, i = 0;
            for(std::map<unsigned int, Eigen::VectorXd>::iterator it_i = knot.begin();
                it_i != knot.end(); it_i ++){
                j = 0; aux = 0.0;
                for(std::map<unsigned int, Eigen::VectorXd>::iterator it_j = knot.begin();
                it_j != knot.end(); it_j ++){
                    aux += iPsi(j,i)*bvp.rbf.Value(it_j->second-center_vec,point-center_vec,c2);
                    j++;
                }
                GN[it_i->first] += -1.0*Y*aux;
                GGN[it_i->first] += pow(Y*aux,2.0);
                i++;
            }
        }
    }

     void  Test_Center_Interpolator(Eigen::VectorXd point, BVP bvp){
        double angle;
        switch(interior){
            case 1:
                angle = Atan0(point);
                break;
            case -1:
                angle = Atan0(point);
                break;
            case -2:
                angle = Atan90(point);
                break;
            case -3:
                angle = Atan180(point);
                break;
            case -4:
                angle = Atan270(point);
                break;
        }
        double aux;
        if(interior == 1){
            for(std::map<unsigned int,double>::iterator it = theta.begin();
                it != theta.end(); it ++){
                aux = 0.0;
                for(unsigned int l = 0;l < N_knots[1]; l++){
                    aux += cos(k[l]*(angle-(it->second)));
                }
                GN[it->first] += (1.0/N_knots[1])*aux;
            }
        } else {
            int j = 0, i = 0;
            for(std::map<unsigned int, Eigen::VectorXd>::iterator it_i = knot.begin();
                it_i != knot.end(); it_i ++){
                j = 0; aux = 0.0;
                for(std::map<unsigned int, Eigen::VectorXd>::iterator it_j = knot.begin();
                it_j != knot.end(); it_j ++){
                    aux += iPsi(j,i)*bvp.rbf.Value(it_j->second,point,c2);
                    j++;
                }
                GN[it_i->first] += aux;
                i++;
            }
        }
    }
    
    void Return_G(unsigned int N_tray, std::vector<double> & G, std::vector<double> & GG, std::vector<int> &  G_j){
        double norm = 1.0/N_tray;
        G.resize(N_knots[0]);
        GG.resize(N_knots[0]);
        G_j.resize(N_knots[0]);
        unsigned int cent = 0;
        for(std::map<unsigned int,double>::iterator it = GN.begin();
            it != GN.end(); it ++){
            G_j[cent] = (int) it->first;
            G[cent] = norm * (it->second);
            GG[cent] = norm * norm * (GGN[it->first]);
            cent ++;
        }
    }
    void Print(void){
        char fname[256];
        sprintf(fname,"Python/center_plot_%u%u.csv",index[0],index[1]);
        FILE *pf;
        pf = fopen(fname,"w");
        fprintf(pf,"x,y,i\n");
        for(std::map<unsigned int, Eigen::VectorXd>::iterator map_it = knot.begin();
            map_it != knot.end();
            map_it ++){
            fprintf(pf,"%f,%f,%u\n",(map_it->second)[0], (map_it->second)[1], map_it->first);
        }
        fclose(pf);
        //char summon_Python[256];
        //sprintf(summon_Python,"python3 center_plot.py");
        //system(summon_Python);
        //getchar();
    }
    
    //Sends the circle information to a worker
    void Send_To_Worker(MPI_Status & status, MPI_Comm & world, unsigned int type){
        if(type == 0){
            MPI_Send(parameters, 3, MPI_DOUBLE, status.MPI_SOURCE, SUPP_CIRCLE_PARAMS, world);
            //printf("Parameters sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(index, 2, MPI_UNSIGNED, status.MPI_SOURCE, SUPP_CIRCLE_INDEX, world);
            //printf("Index is sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(N_knots, 2, MPI_UNSIGNED, status.MPI_SOURCE, SUPP_CIRCLE_N_KNOTS, world);
            //printf("N_knots is sent from server to %d\n", status.MPI_SOURCE);
            double *aux_double  = new double[N_knots[0]*3 + 1];
            unsigned int *aux_uint = new unsigned int[N_knots[0]];
            int *aux_int =  new int[N_knots[1]+ 1];
            int cent = 0;
            for(std::map<unsigned int,Eigen::VectorXd>::iterator it = knot.begin();
                it != knot.end();
                it ++){
                aux_uint[cent] = it -> first;
                aux_double[cent] = (it->second)[0];
                aux_double[N_knots[0] + cent] = (it->second)[1];
                aux_double[2*N_knots[0] + cent] = theta[it->first];
                cent ++;
            }
            for(unsigned int i = 0; i < N_knots[1]; i++){
                aux_int[i] = k[i];
            }
            aux_int[N_knots[1]] = interior;
            aux_double[N_knots[0]*3] = c2;
            MPI_Send(aux_double, N_knots[0]*3 + 1, MPI_DOUBLE, status.MPI_SOURCE, SUPP_CIRCLE_DOUBLES, world);
            //printf("Doubles are sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(aux_uint, N_knots[0], MPI_UNSIGNED, status.MPI_SOURCE, SUPP_CIRCLE_UINTS, world);
            //printf("Uints are sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(aux_int, N_knots[1] + 1, MPI_INT, status.MPI_SOURCE, SUPP_CIRCLE_INTS, world);
            //printf("Ints sent from server to %d\n", status.MPI_SOURCE);
            delete aux_uint; delete aux_double; delete aux_int;
        } else {
            MPI_Send(parameters, 3, MPI_DOUBLE, status.MPI_SOURCE, INT_CIRCLE_PARAMS, world);
            //printf("Parameters sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(index, 2, MPI_UNSIGNED, status.MPI_SOURCE, INT_CIRCLE_INDEX, world);
            //printf("Index is sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(N_knots, 2, MPI_UNSIGNED, status.MPI_SOURCE, INT_CIRCLE_N_KNOTS, world);
            //printf("N_knots is sent from server to %d\n", status.MPI_SOURCE);
            double *aux_double  = new double[N_knots[0]*3 + 1];
            unsigned int *aux_uint = new unsigned int[N_knots[0]];
            int *aux_int =  new int[N_knots[1]+ 1];
            int cent = 0;
            for(std::map<unsigned int,Eigen::VectorXd>::iterator it = knot.begin();
                it != knot.end();
                it ++){
                aux_uint[cent] = it -> first;
                aux_double[cent] = (it->second)[0];
                aux_double[N_knots[0] + cent] = (it->second)[1];
                aux_double[2*N_knots[0] + cent] = theta[it->first];
                cent ++;
            }
            for(unsigned int i = 0; i < N_knots[1]; i++){
                aux_int[i] = k[i];
            }
            aux_double[N_knots[0]*3] = c2;
            aux_int[N_knots[1]] = interior;
            MPI_Send(aux_double, N_knots[0]*3 + 1, MPI_DOUBLE, status.MPI_SOURCE, INT_CIRCLE_DOUBLES, world);
            //printf("Doubles are sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(aux_uint, N_knots[0], MPI_UNSIGNED, status.MPI_SOURCE, INT_CIRCLE_UINTS, world);
            //printf("Uints are sent from server to %d\n", status.MPI_SOURCE);
            MPI_Send(aux_int, N_knots[1] + 1, MPI_INT, status.MPI_SOURCE, INT_CIRCLE_INTS, world);
            //printf("Ints sent from server to %d\n", status.MPI_SOURCE);
            delete aux_uint; delete aux_double; delete aux_int;
        }
        
    }
    //Receives circle information from the server
    void Recieve_From_Server( int server, MPI_Comm & world, unsigned int type, BVP bvp){
        int myid;
        MPI_Comm_rank(world, &myid);
        knot.clear();
        theta.clear();
        k.clear();
        GN.clear();
        GGN.clear();
        if(type == 0){
            MPI_Recv(parameters, 3, MPI_DOUBLE, server, SUPP_CIRCLE_PARAMS, world, MPI_STATUS_IGNORE);
            //printf("Parameters received from server in %d\n", myid);
            //printf("%f %f %f\n",parameters[0],parameters[1], parameters[2]);
            MPI_Recv(index, 2, MPI_UNSIGNED, server, SUPP_CIRCLE_INDEX, world, MPI_STATUS_IGNORE);
            //printf("Index is received from server in %d\n", myid);
            MPI_Recv(N_knots, 2, MPI_UNSIGNED, server, SUPP_CIRCLE_N_KNOTS, world, MPI_STATUS_IGNORE);
            //printf("N_knots is received from server in %d\n", myid);
            //printf("%u %u\n",N_knots[0],N_knots[1]);
            double *aux_double;
            aux_double = new double[N_knots[0]*3 + 1];
            unsigned int *aux_uint;
            aux_uint = new unsigned int[N_knots[0]];
            int *aux_int;
            aux_int = new int[N_knots[1] + 1];
            MPI_Recv(aux_double, N_knots[0]*3 + 1, MPI_DOUBLE, server, SUPP_CIRCLE_DOUBLES, world,MPI_STATUS_IGNORE);
            //printf("Doubles are received from server in %d\n", myid);
            MPI_Recv(aux_uint, N_knots[0], MPI_UNSIGNED, server, SUPP_CIRCLE_UINTS, world,MPI_STATUS_IGNORE);
            //printf("Uints are received from server in %d\n", myid);
            MPI_Recv(aux_int, N_knots[1] + 1, MPI_INT, server, SUPP_CIRCLE_INTS, world, MPI_STATUS_IGNORE);
            //printf("Ints  are received from server in %d\n", myid);
            Eigen::VectorXd aux_vec;
            aux_vec.resize(2);
            k.resize(N_knots[1]);
            for(unsigned int i = 0;
                i < N_knots[0];
                i++){
                aux_vec[0] = aux_double[i];
                aux_vec[1] = aux_double[N_knots[0] + i];
                knot[aux_uint[i]] = aux_vec;
                theta[aux_uint[i]] = aux_double[2*N_knots[0] + i];
                GN[aux_uint[i]] = 0.0;
                GGN[aux_uint[i]] = 0.0;
            }
            for(unsigned int i = 0; i < N_knots[1]; i++){
                k[i] = aux_int[i];
            }
            center_vec.resize(2);
            center_vec[0] = parameters[1];
            center_vec[1] = parameters[2];
            c2 = aux_double[N_knots[0]*3];
            int i,j;
            Eigen::MatrixXd Psi, I;
            Psi.resize(N_knots[0],N_knots[0]);
            I = Psi;
            I.setIdentity();
            i = 0;
            for(std::map<unsigned int, Eigen::VectorXd>::iterator it_i = knot.begin();
                it_i != knot.end(); it_i ++){
                    j = 0;
                    for(std::map<unsigned int, Eigen::VectorXd>::iterator it_j = knot.begin();
                    it_j != knot.end(); it_j ++){
                        Psi(i,j) = bvp.rbf.Value(it_i->second-center_vec,it_j->second-center_vec,c2);
                        if(isnan(Psi(i,j))) printf("THERE IS A NAN IN PSI\n");
                        j++;
                    }
                    i++;
            }
            Eigen::FullPivLU<Eigen::MatrixXd> lu(Psi);
            if(lu.isInvertible()){
                iPsi = lu.inverse();
            }else{
                printf("FATAL ERROR: Psi Matrix is not invertible.");
                iPsi = Psi*0.0;
            }
            interior = aux_int[N_knots[1]];

            delete aux_double; delete aux_uint; delete aux_int;
        }else{
            MPI_Recv(parameters, 3, MPI_DOUBLE, server, INT_CIRCLE_PARAMS, world, MPI_STATUS_IGNORE);
            //printf("Parameters received from server in %d\n", myid);
            MPI_Recv(index, 2, MPI_UNSIGNED, server, INT_CIRCLE_INDEX, world, MPI_STATUS_IGNORE);
            //printf("Index is received from server in %d\n", myid);
            //printf("Int index %u %u\n",index[0],index[1]);
            MPI_Recv(N_knots, 2, MPI_UNSIGNED, server, INT_CIRCLE_N_KNOTS, world, MPI_STATUS_IGNORE);
            //printf("N_knots is received from server in %d\n", myid);
            double *aux_double;
            aux_double = new double[N_knots[0]*3 + 1];
            unsigned int *aux_uint;
            aux_uint = new unsigned int[N_knots[0]];
            int *aux_int;
            aux_int = new int[N_knots[1] + 1];
            MPI_Recv(aux_double, N_knots[0]*3 + 1, MPI_DOUBLE, server, INT_CIRCLE_DOUBLES, world,MPI_STATUS_IGNORE);
            //printf("Doubles are received from server in %d\n", myid);
            MPI_Recv(aux_uint, N_knots[0], MPI_UNSIGNED, server, INT_CIRCLE_UINTS, world,MPI_STATUS_IGNORE);
            //printf("Uints are received from server in %d\n", myid);
            MPI_Recv(aux_int, N_knots[1] + 1, MPI_INT, server, INT_CIRCLE_INTS, world, MPI_STATUS_IGNORE);
            //printf("Ints  are received from server in %d\n", myid);
            Eigen::VectorXd aux_vec;
            aux_vec.resize(2);
            k.resize(N_knots[1]);
            for(unsigned int i = 0;
                i < N_knots[0];
                i++){
                //printf("%u\n", aux_uint[i]);
                aux_vec[0] = aux_double[i];
                aux_vec[1] = aux_double[N_knots[0] + i];
                knot[aux_uint[i]] = aux_vec;
                theta[aux_uint[i]] = aux_double[2*N_knots[0] + i];
                GN[aux_uint[i]] = 0.0;
                GGN[aux_uint[i]] = 0.0;
            }
            for(unsigned int i = 0; i < N_knots[1]; i++){
                k[i] = aux_int[i];
            }
            center_vec.resize(2);
            center_vec[0] = parameters[1];
            center_vec[1] = parameters[2];
            c2 = aux_double[N_knots[0]*3];
            int i,j;
            Eigen::MatrixXd Psi, I;
            Psi.resize(N_knots[0],N_knots[0]);
            I = Psi;
            I.setIdentity();
            i = 0;
            for(std::map<unsigned int, Eigen::VectorXd>::iterator it_i = knot.begin();
                it_i != knot.end(); it_i ++){
                    j = 0;
                    for(std::map<unsigned int, Eigen::VectorXd>::iterator it_j = knot.begin();
                    it_j != knot.end(); it_j ++){
                        Psi(i,j) = bvp.rbf.Value(it_i->second-center_vec,it_j->second-center_vec,c2);
                        if(isnan(Psi(i,j))) printf("THERE IS A NAN IN PSI\n");
                        j++;
                    }
                    i++;
            }
            Eigen::FullPivLU<Eigen::MatrixXd> lu(Psi);
            if(lu.isInvertible()){
                iPsi = lu.inverse();
            }else{
                printf("FATAL ERROR: Psi Matrix is not invertible.");
                iPsi = Psi*0.0;
            }
            interior = aux_int[N_knots[1]];
            delete aux_double; delete aux_uint; delete aux_int;
        } 
    }
};
#endif
