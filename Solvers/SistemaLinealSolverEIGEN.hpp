#ifndef SISTEMALINEALSOLVEREIGEN
#define SISTEMALINEALSOLVEREIGEN
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "../Mesh/SistemaLineal.hpp"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
typedef Eigen::Triplet<double,int> Triplet;
class SistemaLinealSolverEIGEN{
    public:
    std::vector<double> solucion, residuo;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> G;
    /*SistemaLinealSolverEIGEN(){
        solucion.clear();
        residuo.clear();
        Eigen::SparseMatrix<double> G;
        G.resize(0,0);
        solver.factorize(G);
    }*/
    void Initialize(SistemaLineal & sistema){
        G.resize(sistema.n_filas,sistema.n_filas);
        G.reserve(sistema.nnz);
        std::vector<Triplet> G_triplet;
        G_triplet.resize(sistema.nnz);
        std::vector<int> G_i,G_j;
        std::vector<double> G_ij;
        sistema.Convierte_COO(G_i,G_j,G_ij);
        for(int i = 0; i < sistema.nnz; i++){
            if(fabs(G_ij[i])>1E-20) G_triplet[i] = Triplet(G_i[i], G_j[i], G_ij[i]);
        }
        G.setFromTriplets(G_triplet.begin(),G_triplet.end());
        G_triplet.clear();
        // Compute the ordering permutation vector from the structural pattern of A
        //solver.analyzePattern(G); 
        // Compute the numerical factorization 
        //solver.factorize(G); 
        //Use the factors to solve the linear system 
    }
    void Solve(SistemaLineal & sistema){
        Eigen::VectorXd B, sol;
        B.resize(sistema.n_filas);
        for(int i = 0; i < sistema.n_filas; i++){
            B[i] = sistema.B[i];
        }
        sol = solver.solve(B);
        solucion.resize(sistema.n_filas);
        residuo.resize(sistema.n_filas);
        for(int i = 0; i < sistema.n_filas; i++){
            solucion[i] = sol[i];
            sistema.u[i] = sol[i];
        }
    };
    void Estima_Numero_Condicionamiento(int N_max, int N_min){
        Eigen::SparseMatrix<double,Eigen::RowMajor> GG;
        std::cout << "Computing G*G'..." << std::endl;
        GG = G*G.transpose();
        Eigen::VectorXd vector;
        vector.resize(GG.rows());
        for(int i = 0; i < vector.rows(); i++){
            vector[i] = (i%2==0) ? 1.0 : -1.0;
        }
        std::cout << "Computing Max singular value..." << std::endl;
        for(int i = 0; i < N_max; i++){
            vector = vector/vector.norm();
            vector = GG*vector;
        }
        std::cout << "Max singular value " << sqrt(vector.norm()) << std::endl;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::IncompleteLUT<double> > BiCGSTAB;
        std::cout << "Computing GG' BiCGSTAB solver" << std::endl;
        solver.analyzePattern(GG); 
        solver.factorize(G); 
        //solver.compute(GG);
        for(int i = 0; i < vector.rows(); i++){
            vector[i] = 1.0;
        }
        std::cout << "Computing Min singular value..." << std::endl;
        for(int i = 0; i < N_min; i++){
            //vector = vector/vector.norm();
            vector = solver.solve(vector/vector.norm());
            std::cout << i << " " << 1/sqrt(vector.norm()) << std::endl;
        }
        std::cout << "Min singular value " << 1/sqrt(vector.norm()) << std::endl;
    }
};
#endif
