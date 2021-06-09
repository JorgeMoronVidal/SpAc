#ifndef SPAC
#define SPAC
#include<eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include "GMSolver.hpp"
#include "BVP.hpp"
#include <algorithm>
#define ASK_FOR_JOB 100
#define REPLY_WORKER 101
#define TAG_Gi 120
#define TAG_Gj 121
#define TAG_Gval 122
#define TAG_Gvar 122
#define TAG_Bi 130
#define TAG_Bval 131
#define TAG_Bvar 122
#define K_EXCESS 1.9f
//#define AVERAGED
typedef Eigen::Triplet<double,int> T;
class SpAc_solver{
    private:
    /*
    MPI-related variables
    */
    MPI_Comm world, workers;
    MPI_Group world_group, worker_group;
    MPI_Status status;
    int ranks[1],work_control[2],numprocs, myid, server, workerid;
    /*Filename of the Flux DEBUG File*/
    char debug_fname[256];
    //Vector of the centers considered
    std::vector<Center> centers;
    std::map<std::vector<Center>::iterator, bool> center_solved;
    //Vector of the Knots and its support [0 index] integration centers from [1] in advance
    std::vector<std::vector<std::vector<Center>::iterator> > knots;
    double domain_parameters[4];
    unsigned int N_knots;
    /*G and B storage vector*/
    std::vector<double> G, var_G, B, var_B;
    std::vector<int>  G_j, G_i, B_i;
    /*Triplet's vector*/
    std::vector<T> T_vec_G, T_vec_B;
    public:
    SpAc_solver(int argc, char *argv[]){
        domain_parameters[0] = domain_parameters[1] = domain_parameters[2] = domain_parameters[3] = 0.0;
        MPI_Configuration(argc, argv);
    }
    void MPI_Configuration(int argc, char *argv[]){
            MPI_Init(&argc, &argv);
            world = MPI_COMM_WORLD;
            MPI_Comm_size(world, &numprocs);
            MPI_Comm_rank(world, &myid);
            server = numprocs-1;
            MPI_Comm_group(world, &world_group);
            ranks[0] = server;
            //We exclude a group in world to create workers
            MPI_Group_excl(world_group, 1, ranks, &worker_group);
            //Workers comunicator is created
            MPI_Comm_create(world, worker_group, &workers);
            //We free both comunicators groups
            MPI_Group_free(&worker_group);
            MPI_Group_free(&world_group);
            sprintf(debug_fname,"Debug/Debug_output_p%d.txt",myid);
            FILE *dfile;
            dfile = fopen(debug_fname,"w");
            fclose(dfile);
        }
    void Init(double global_parameters[4], unsigned int N_centers, double superposition_coefficent, 
            unsigned int N_knots_per_circle, BVP bvp, double fac){
            double d,R, disp;
            Center aux_center;
            unsigned int knot_centinel = 0;
            domain_parameters[0] = global_parameters[0];
            domain_parameters[1] = global_parameters[1];
            domain_parameters[2] = global_parameters[2];
            domain_parameters[3] = global_parameters[3];
            disp = (domain_parameters[2] - domain_parameters[0])*superposition_coefficent/(2*superposition_coefficent + 2*K_EXCESS*(N_centers-1));
            d = (domain_parameters[2] - domain_parameters[0] - 2*disp)/(N_centers -1);
            if(superposition_coefficent > 2.0){
                printf("FATAL WARNING: Superposition coefficent has to be equal or lower to 2.0\n");
            }
            if(superposition_coefficent < 1.0){
                printf("FATAL WARNING: Superposition coefficent has to be equal or higher than 1.0\n");
            }
            R = superposition_coefficent*d/sqrt(2.0);
            if(myid == server){
            system("rm Python/*");
            system("rm Debug/*");
            system("rm Matlab_buffer/*");
            Eigen::VectorXd aux_vec;
            aux_vec.resize(2);
            //Initialization of all the centers
            for(unsigned int h_index = 0; h_index <= N_centers-1; h_index++){
                for(unsigned int v_index = 0; v_index <= N_centers-1; v_index++){
                    aux_vec[0] = global_parameters[0] + disp + h_index*d; aux_vec[1] = global_parameters[1] + disp + v_index*d;
                    knot_centinel = aux_center.Init_Square(R,bvp, domain_parameters,fac,aux_vec,h_index,v_index,N_knots_per_circle,knot_centinel,N_centers-1,N_centers-1);
                    centers.push_back(aux_center);
                    //printf("%u %u\n",h_index,v_index);
                    aux_center.Print();
                }
            }
            N_knots = knot_centinel;
            knots.resize(N_knots);
            char fname[256];
            FILE *file;
            //Each interior center has a file dedicated to it
            for(std::vector<Center>::iterator center_it = centers.begin();
                center_it != centers.end();
                center_it ++){
                center_solved[center_it] = false;
                if((*center_it).Is_interior()){
                    sprintf(fname,"Matlab_buffer/patch_%u%u.csv",(*center_it).index[0],(*center_it).index[1]);
                    file = fopen(fname,"w");
                    fclose(file);
                }
            }
            double d_centinel;
            //The center to which each knot belongs is identified and stored in the first element of knots
            for(unsigned int knot_index = 0; knot_index < N_knots; knot_index ++){
                for(std::vector<Center>::iterator center_it = centers.begin();
                center_it != centers.end();
                center_it ++){
                    if((*center_it).Get_knot_position(knot_index,aux_vec)){
                        knots[knot_index].push_back(center_it);
                        printf("Knot %u is in center [%u %u]\n",knot_index,(*center_it).index[0],(*center_it).index[1]);
                    }
                }
                //The center or centers in whic each knot is integrated are identified and stored the the latter elements 
                #ifdef AVERAGED
                for(std::vector<Center>::iterator center_it = centers.begin();
                center_it != centers.end();
                center_it ++){
                    if((*center_it).Is_inside(aux_vec)){
                        knots[knot_index].push_back(center_it);
                        if((*center_it).Is_interior()){
                            sprintf(fname,"Matlab_buffer/patch_%u%u.csv",(*center_it).index[0],(*center_it).index[1]);
                            file = fopen(fname,"a");
                            fprintf(file,"%u,%f,%f\n",knot_index,aux_vec[0],aux_vec[1]);
                            fclose(file);
                        }
                        printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*center_it).index[0],(*center_it).index[1]);
                    }
                }
                #endif
                #ifndef AVERAGED
                knots[knot_index].push_back(centers.begin());
                d_centinel = 0.0;
                for(std::vector<Center>::iterator center_it = centers.begin();
                center_it != centers.end();
                center_it ++){
                    if((*center_it).Is_inside(aux_vec)){
                        if((*center_it).Boundary_Distance(aux_vec) < d_centinel){
                            d_centinel = (*center_it).Boundary_Distance(aux_vec);
                            knots[knot_index][1] = center_it;
                        }
                    }
                }
                if((*knots[knot_index][1]).Is_interior()){
                    sprintf(fname,"Matlab_buffer/patch_%u%u.csv",(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                    file = fopen(fname,"a");
                    fprintf(file,"%u,%f,%f\n",knot_index,aux_vec[0],aux_vec[1]);
                    fclose(file);
                }
                printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                #endif
            }
            system("python3 center_plot.py");
        }
    }
    void Send_G_B(void){
    work_control[0] = (int) G.size();
    work_control[1] = (int) B.size();
    MPI_Send(work_control, 2, MPI_INT, server, ASK_FOR_JOB, world);
    //G_i is sent
    int *G_i_aux;
    G_i_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_i_aux[i] = G_i[i];
    G_i.resize(0);
    MPI_Send(G_i_aux, work_control[0], MPI_INT, server, TAG_Gi, world);
    delete G_i_aux;
    //G_j is sent
    int *G_j_aux;
    G_j_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_j_aux[i] = G_j[i];
    G_j.resize(0);
    MPI_Send(G_j_aux, work_control[0], MPI_INT, server, TAG_Gj, world);
    delete G_j_aux;
    //G_val is sent
    double *G_aux;
    G_aux = new double[work_control[0]];
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G[i];
    }
    G.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gval, world);
    delete G_aux;
    //B_i is sent
    int *B_i_aux;
    B_i_aux = new int[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_i_aux[i] = B_i[i];
    B_i.resize(0);
    MPI_Send(B_i_aux, work_control[1], MPI_INT, server, TAG_Bi, world);
    delete B_i_aux;
    //B is sent
    double *B_aux;
    B_aux = new double[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B[i];
    B.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bval, world);
    delete B_aux;
    FILE *dfile;
    dfile = fopen(debug_fname,"a");
    fprintf(dfile,"Processor %d ended sending its contribution to G and B\n", myid);
    fclose(dfile);
    MPI_Comm_free(& workers);
}
    void Receive_G_B(void){
    //It ends worker process solving
    MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
    ASK_FOR_JOB, world, &status);
    work_control[0] = 0;
    work_control[1] = 0;
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
    //The server accumulates the values of G and B
    MPI_Recv(work_control, 2, MPI_INT, status.MPI_SOURCE, ASK_FOR_JOB, world, &status);
    G.resize(work_control[0]);
    G_i.resize(work_control[0]);
    G_j.resize(work_control[0]);
    //G_i is received
    int *G_i_aux, *G_j_aux;
    double *G_aux;
    G_i_aux = new int[work_control[0]];
    MPI_Recv(G_i_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gi, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_i[i] = G_i_aux[i];
    delete G_i_aux;
    //G_j is received
    G_j_aux = new int[work_control[0]];
    MPI_Recv(G_j_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gj, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_j[i] = G_j_aux[i];
    delete G_j_aux;
    //G is received
    G_aux = new double[work_control[0]];
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gval, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    delete G_aux;
    //Triplets vector is fullfilled
    for(int j = 0; j < work_control[0]; j++){
        T_vec_G.push_back(T(G_i[j], G_j[j], G[j]));
    }
    G.resize(0);
    G_i.resize(0);
    G_j.resize(0);
    B.resize(work_control[1]);
    B_i.resize(work_control[1]);
    double *B_aux;
    int *B_i_aux;
    //B_i is received
    B_i_aux = new int[work_control[1]];
    MPI_Recv(B_i_aux, work_control[1], MPI_INT, status.MPI_SOURCE, TAG_Bi, world, &status);
    for(int i = 0; i< work_control[1]; i++) B_i[i] = B_i_aux[i];
    delete B_i_aux;
    //B is received
    B_aux = new double[work_control[1]];
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bval, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    delete B_aux;
    for(int j = 0; j < work_control[1]; j++){
        T_vec_B.push_back(T(B_i[j], 0, B[j]));
    }
}
    void Compute_Solution(BVP bvp){
        Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
        //Eigen::SparseMatrix<double> I_sparse;
        unsigned int size = N_knots;
        B_sparse.resize(size, 1);
        B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
        T_vec_B.resize(0);
        //I_sparse.resize(size, size);
        //I_sparse.setIdentity();
        G_sparse.resize(size, size);
        G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
        T_vec_G.resize(0);
        //G_sparse+=I_sparse;
        for(unsigned int i = 0; i < knots.size(); i ++) G_sparse.coeffRef(i,i) =  std::max((double)(knots[i].size()-1),1.0);
        I_sparse.resize(0,0);
        FILE *ofile;
        ofile = fopen("Debug/G.csv","w");
        fprintf(ofile, "G_i,G_j,G_ij\n");
        for (int k=0; k<G_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld,%ld,%.8e\n",it.row(),it.col(),it.value());
            }
        fclose(ofile);
        ofile = fopen("Debug/B.csv","w");
        fprintf(ofile, "B_i,B_j,B_ij\n");
        for (int k=0; k<B_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld,%ld,%.8e\n",it.row(),it.col(),it.value());
            }
    fclose(ofile);
    Eigen::VectorXd ud,Bd;
    Bd = Eigen::VectorXd(B_sparse);
    B_sparse.resize(0,0);
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_LU;
    solver_LU.compute(G_sparse);
    solver_LU.factorize(G_sparse);
    ud = solver_LU.solve(Bd);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_IT;
    solver_IT.compute(G_sparse);
    ud = solver_IT.solveWithGuess(Bd,ud);
    ofile = fopen("Output/solution.csv","w");
    fprintf(ofile,"Knot_index,x,y,sol_analytic,sol_PDDS,err,rerr\n");
    double sol, err, rerr;
    Eigen::VectorXd Knot_Position;
    for(unsigned int i = 0; i<ud.size(); i++){
        (*knots[i][0]).Get_knot_position(i,Knot_Position);
        sol = bvp.u.Value(Knot_Position,INFINITY);
        err = ud(i) - sol;
        rerr = fabs(err)/sol;
        fprintf(ofile,"%u,%1.10f,%1.6f,%1.6f,%1.10f,%1.10f,%1.10f\n",i,Knot_Position[0],Knot_Position[1],
        sol,ud(i),err,rerr);
    }
    fclose(ofile);
}
    void Write_G_B(char file_G_i[256], char file_G_j[256], char file_G[256],
                   char file_B_i[256], char file_B[256]){
        Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
        //Eigen::SparseMatrix<double> I_sparse;
        unsigned int size = N_knots;
        B_sparse.resize(size, 1);
        B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
        T_vec_B.resize(0);
        //I_sparse.resize(size, size);
        //I_sparse.setIdentity();
        G_sparse.resize(size, size);
        G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
        T_vec_G.resize(0);
        //G_sparse+=I_sparse;
        //for(unsigned int i = 0; i < knots.size(); i ++) G_sparse.coeffRef(i,i) =  std::max((double)(knots[i].size()-1),1.0);
        I_sparse.resize(0,0);
        FILE *ofile;
        ofile = fopen(file_G_i,"w");
        for (int k=0; k<G_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld\n",it.row());
            }
        fclose(ofile);
        ofile = fopen(file_G_j,"w");
        for (int k=0; k<G_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld\n",it.col());
            }
        fclose(ofile);
        ofile = fopen(file_G,"w");
        for (int k=0; k<G_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
            {
                fprintf(ofile,"%.8e\n",it.value());
            }
        fclose(ofile);
        ofile = fopen(file_B_i,"w");
        for (int k=0; k<B_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld\n",it.row());
            }
        fclose(ofile);
        ofile = fopen(file_B,"w");
        for (int k=0; k<B_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
            {
                fprintf(ofile,"%.8e\n",it.value());
            }
        fclose(ofile);
    }
    void Read_G_B_And_Solve(BVP bvp, char G_i_MC[256], char G_j_MC[256], char G_MC[256], char B_i_MC[256],
        char B_MC[256], char G_i_Det[256], char G_j_Det[256], char G_Det[256], char B_i_Det[256], char B_Det[256]){
        if(myid == server){
        G.resize(0);
        G_i.resize(0);
        G_j.resize(0);
        B.resize(0);
        B_i.resize(0);
        FILE *data;
        data = fopen(G_i_MC,"r");
        char *buffer = NULL;
        size_t line_buf_size;
        ssize_t line_size;
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_i.push_back(atoi(buffer));
            printf("G_i_MC %d \n",atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(G_i_Det,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_i.push_back(atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(G_j_MC,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_j.push_back(atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(G_j_Det,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_j.push_back(atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(G_MC,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G.push_back(atof(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(G_Det,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G.push_back(atof(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        //Triplets vector is fullfilled
        for(unsigned int j = 0; j <G.size(); j++){
            T_vec_G.push_back(T(G_i[j], G_j[j], G[j]));
        }
        G.resize(0);
        G_i.resize(0);
        G_j.resize(0);
        data = fopen(B_i_MC,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B_i.push_back(atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(B_i_Det,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B_i.push_back(atoi(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(B_MC,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B.push_back(atof(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        data = fopen(B_Det,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B.push_back(atof(buffer));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        for(unsigned int j = 0; j <B_i.size(); j++){
            T_vec_B.push_back(T(B_i[j], 0, B[j]));
        }
        B.resize(0);
        B_i.resize(0);
        Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
        //Eigen::SparseMatrix<double> I_sparse;
        unsigned int size = N_knots;
        B_sparse.resize(size, 1);
        B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
        T_vec_B.resize(0);
        //I_sparse.resize(size, size);
        //I_sparse.setIdentity();
        G_sparse.resize(size, size);
        G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
        T_vec_G.resize(0);
        //G_sparse+=I_sparse;
        for(unsigned int i = 0; i < knots.size(); i ++) G_sparse.coeffRef(i,i) =  std::max((double)(knots[i].size()-1),1.0);
        I_sparse.resize(0,0);
        FILE *ofile;
        ofile = fopen("Debug/G.csv","w");
        fprintf(ofile, "G_i,G_j,G_ij\n");
        for (int k=0; k<G_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld,%ld,%.8e\n",it.row(),it.col(),it.value());
            }
        fclose(ofile);
        ofile = fopen("Debug/B.csv","w");
        fprintf(ofile, "B_i,B_j,B_ij\n");
        for (int k=0; k<B_sparse.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
            {
                fprintf(ofile,"%ld,%ld,%.8e\n",it.row(),it.col(),it.value());
            }
        fclose(ofile);
        Eigen::VectorXd ud,Bd;
        Bd = Eigen::VectorXd(B_sparse);
        B_sparse.resize(0,0);
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_LU;
        solver_LU.compute(G_sparse);
        solver_LU.factorize(G_sparse);
        ud = solver_LU.solve(Bd);
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_IT;
        solver_IT.compute(G_sparse);
        ud = solver_IT.solveWithGuess(Bd,ud);
        ofile = fopen("Output/solution.csv","w");
        fprintf(ofile,"Knot_index,x,y,sol_analytic,sol_PDDS,err,rerr\n");
        double sol, err, rerr;
        Eigen::VectorXd Knot_Position;
        for(unsigned int i = 0; i<ud.size(); i++){
            (*knots[i][0]).Get_knot_position(i,Knot_Position);
            sol = bvp.u.Value(Knot_Position,INFINITY);
            err = ud(i) - sol;
            rerr = fabs(err)/sol;
            fprintf(ofile,"%u,%1.10f,%1.6f,%1.6f,%1.10f,%1.10f,%1.10f\n",i,Knot_Position[0],Knot_Position[1],
            sol,ud(i),err,rerr);
        }
        fclose(ofile);
        }
    }
    void Test_Interpolator(BVP bvp){
        FILE *pFile;
        Eigen::VectorXd point, aux_vec;
        double angle, value;
        if(myid == server){
            pFile = fopen ("Debug/Interpolator.csv","w");
            fprintf(pFile, "index,center,x,y,u_exact,u_interp\n");
            fclose(pFile);
            for(unsigned int knot_index = 0; knot_index < knots.size(); knot_index ++){
                MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_FOR_JOB, world, &status);
                work_control[0] = 1;
                work_control[1] = knot_index;
                MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
                //printf("Knot %u is in center [%u %u]\n",knot_index,(*knots[knot_index][0]).index[0],(*knots[knot_index][0]).index[1]);
                //printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                (*knots[knot_index][0]).Send_To_Worker(status, world,0);
                (*knots[knot_index][1]).Send_To_Worker(status, world,1);
            }
            for(int worker_index = 0; worker_index < server; worker_index ++){
                MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                ASK_FOR_JOB, world, &status);
                work_control[0] = 0;
                work_control[1] = 0;
                MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
            }
        } else {
            double start = MPI_Wtime();;
            work_control[0] = 1;
            Center supp_center, int_center;
            while(work_control[0] == 1){
                MPI_Send(work_control, 2, MPI_INT, server, ASK_FOR_JOB, world);
                MPI_Recv(work_control, 2, MPI_INT, server, REPLY_WORKER, world, &status);
                if(work_control[0] == 1){
                    supp_center.Recieve_From_Server(server,world,0,bvp);
                    int_center.Recieve_From_Server(server,world,1,bvp);
                    angle = supp_center.theta[work_control[1]] + (std::rand()/RAND_MAX -0.5)*M_PI/48;
                    point.resize(2);
                    point[0] = supp_center.center_vec[0] + supp_center.Get_radius()*cos(angle);
                    point[1] = supp_center.center_vec[1] + supp_center.Get_radius()*sin(angle);
                    supp_center.Test_Center_Interpolator(point,bvp);
                    value = 0.0;
                    for(std::map<unsigned int, double>::iterator it_G = supp_center.GN.begin();
                        it_G != supp_center.GN.end(); it_G ++){
                            aux_vec.resize(2);
                            supp_center.Get_knot_position(it_G->first,aux_vec);
                            value += it_G->second * bvp.u.Value(aux_vec,0.0);
                    }
                    pFile = fopen ("Debug/Interpolator.csv","a");
                    fprintf(pFile, "%u,[%u %u],%f,%f,%f,%f,%e\n",work_control[1],supp_center.index[0],supp_center.index[1],
                    point[0], point[1], bvp.u.Value(point, 0.0), value, bvp.u.Value(point, 0.0) - value);
                    fclose(pFile);
                }
            }
            printf("Process %d is done interpolationg\n", myid);
            double end = MPI_Wtime();
            printf("Process %d used %f  hours \n", myid, (end-start)/3600);
        }
    }; 
    void Solve_ExKnots_MC(BVP bvp, double discretization, unsigned int N_tray){
        double start = MPI_Wtime();
        FILE *pFile;
        if(myid == server){
            pFile = fopen ("Debug/Node_debug.csv","w");
            fprintf(pFile, "index,var_B,APL,time\n");
            fclose(pFile);
            pFile = fopen("Debug/G_var.csv","w");
            if (pFile == NULL) perror("Failed: ");
            fprintf(pFile,"G_i,G_j,var_G_ij\n");
            fclose(pFile);
            for(unsigned int knot_index = 0; knot_index < knots.size(); knot_index ++){
                for(unsigned int index_i_center = 1; index_i_center < knots[knot_index].size(); 
                    index_i_center ++){
                        if(!(*knots[knot_index][index_i_center]).Is_interior()){
                            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_FOR_JOB, world, &status);
                            work_control[0] = 1;
                            work_control[1] = knot_index;
                            MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
                            //printf("Knot %u is in center [%u %u]\n",knot_index,(*knots[knot_index][0]).index[0],(*knots[knot_index][0]).index[1]);
                            //printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                            (*knots[knot_index][0]).Send_To_Worker(status, world,0);
                            (*knots[knot_index][index_i_center]).Send_To_Worker(status, world,1);
                    }
                }
            }
            for(int worker_index = 0; worker_index < server; worker_index ++){
                Receive_G_B();
            }
            char file_G_i[256],file_G_j[256],file_G[256],file_B_i[256], file_B[256];
            sprintf(file_G_i,"Output/G_i_MC.txt");
            sprintf(file_G_j,"Output/G_j_MC.txt");
            sprintf(file_G,"Output/G_MC.txt");
            sprintf(file_B_i,"Output/B_i_MC.txt");
            sprintf(file_B,"Output/B_MC.txt");
            Write_G_B(file_G_i, file_G_j, file_G, file_B_i, file_B);
            pFile = fopen(debug_fname,"a");
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fprintf(pFile,"Process %d ended its work.\n",myid);
            fclose(pFile);
        } else {
            //Start and end time for the node
            double knot_start, knot_end;
            Eigen::VectorXd X0;
            //G and B storage vectors for each node
            double B_temp,B_var_temp;
            std::vector<double> Bv_temp, G_temp, var_G_temp, bparams_vec(domain_parameters, domain_parameters + 4 * sizeof(domain_parameters[0]));
            std::vector<int>  G_i_temp, G_j_temp, B_i_temp;
            //Stencil
            GMSolver solver(bvp ,bparams_vec, discretization, (unsigned int) myid +1);
            work_control[0] = 1;
            Center supp_center, int_center;
            while(work_control[0] == 1){
                MPI_Send(work_control, 2, MPI_INT, server, ASK_FOR_JOB, world);
                MPI_Recv(work_control, 2, MPI_INT, server, REPLY_WORKER, world, &status);
                if(work_control[0] == 1){
                    knot_start = MPI_Wtime();
                    supp_center.Recieve_From_Server(server,world,0,bvp);
                    int_center.Recieve_From_Server(server,world,1,bvp);
                    G_j_temp.resize(0);
                    G_temp.resize(0);
                    supp_center.Get_knot_position(work_control[1],X0);
                    solver.Solve(supp_center,int_center,work_control[1],N_tray,discretization,
                    G_j_temp,G_temp,var_G_temp,B_temp,B_var_temp);
                    B.push_back(B_temp);
                    //B.push_back(bvp.u.Value(X0,0.0));
                    B_i.push_back((int)work_control[1]);
                    var_B.push_back(B_var_temp);
                    for(unsigned int k = 0;k < G_temp.size(); k ++){
                        G_i.push_back((int)work_control[1]);
                        G_j.push_back(G_j_temp[k]);
                        G.push_back(G_temp[k]);
                    }    
                    knot_end = MPI_Wtime();
                    pFile = fopen("Debug/Node_debug.csv","a");
                    fprintf(pFile,"%d,%f,%f,%f\n",work_control[1],
                    solver.var_B,solver.APL,(knot_end-knot_start));
                    fclose(pFile);
                    pFile = fopen("Debug/G_var.csv","a");
                    for(unsigned int i = 0; i < solver.var_G.size(); i++){
                        fprintf(pFile,"%d,%d,%f\n",work_control[1],
                        G_j_temp[i],solver.var_G[i]);
                    }
                    fclose(pFile);
                }
            }
            pFile = fopen(debug_fname,"a");
            fprintf(pFile,"Process %d is done solving nodes \n", myid);
            Send_G_B();
            fprintf(pFile,"Process %d is done sending G \n", myid);
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fclose(pFile);
        }
        printf("Process %d ended Solve_Interfaces_MC\n",myid);
    }
    void Deterministic_Solve(Center supp_center, Center int_center,
        std::vector<int> & G_i_temp, std::vector<int> & G_j_temp, std::vector<double> & G_temp, 
        std::vector<double> & B_temp, std::vector<int> & B_i_temp){
        B_temp.clear();
        B_i_temp.clear();
        G_j_temp.clear();
        G_i_temp.clear();
        G_temp.clear();
        char fname[256];
        sprintf(fname,"Matlab_buffer/parameters_%d.csv",myid);
        FILE *data;
        data = fopen(fname,"w");
        fprintf(data,"%u,%u,%f,%f,%f\n",int_center.index[0],int_center.index[1],
        int_center.Get_center_position()[0],int_center.Get_center_position()[1],
        int_center.Get_radius());
        fclose(data);
        sprintf(fname,"Matlab_buffer/points_%d.csv",myid);
        data = fopen(fname,"w");
        for(std::map<unsigned int, double>::iterator it = int_center.theta.begin();
            it != int_center.theta.end(); it ++){
            fprintf(data,"%u,%f\n",it->first, it->second);
        }
        fclose(data);
        char system_order[512];
        sprintf(system_order,"/usr/local/MATLAB/R2020a/bin/matlab -batch \"cd(\'/home/jorge/Dropbox/DOC/SPDDS\'); Solve_Deterministic(%d);\"",myid);
        system(system_order);
        sprintf(fname,"Matlab_buffer/B_%d.txt",myid);
        //sprintf(fname,"Matlab_buffer/B_0.txt");
        data = fopen(fname,"r");
        char *buffer = NULL;
        size_t line_buf_size;
        ssize_t line_size;
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B_temp.push_back(atof(buffer));
            //printf("G = %f\n", *(--G.end()));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        sprintf(fname,"Matlab_buffer/Bi_%d.txt",myid);
        //sprintf(fname,"Matlab_buffer/Bi_0.txt");
        data = fopen(fname,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            B_i_temp.push_back(atoi(buffer));
            //printf("G = %f\n", *(--G.end()));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        sprintf(fname,"Matlab_buffer/Gi_%d.txt",myid);
        //sprintf(fname,"Matlab_buffer/Gi_0.txt");
        data = fopen(fname,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_i_temp.push_back(atoi(buffer));
            //printf("G = %f\n", *(--G.end()));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        sprintf(fname,"Matlab_buffer/Gj_%d.txt",myid);
        //sprintf(fname,"Matlab_buffer/Gj_0.txt");
        data = fopen(fname,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_j_temp.push_back(atoi(buffer));
            //printf("G = %f\n", *(--G.end()));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
        sprintf(fname,"Matlab_buffer/G_%d.txt",myid);
        //sprintf(fname,"Matlab_buffer/G_0.txt");
        data = fopen(fname,"r");
        line_size = getline(&buffer, &line_buf_size, data);
        while (line_size >= 0){
            G_temp.push_back(atof(buffer));
            //printf("G = %f\n", *(--G.end()));
            line_size = getline(&buffer, &line_buf_size, data);
        }
        fclose(data);
    };
    void Solve_IntKnots_Det(BVP bvp){
        double start = MPI_Wtime();
        FILE *pFile;
        if(myid == server){
            pFile = fopen ("Debug/Node_debug.csv","w");
            fprintf(pFile, "index,var_B,APL,time\n");
            fclose(pFile);
            pFile = fopen("Debug/G_var.csv","w");
            if (pFile == NULL) perror("Failed: ");
            fprintf(pFile,"G_i,G_j,var_G_ij\n");
            fclose(pFile);
            for(unsigned int knot_index = 0; knot_index < knots.size(); knot_index ++){
                for(unsigned int index_i_center = 1; index_i_center < knots[knot_index].size(); 
                    index_i_center ++){
                        if(!center_solved[knots[knot_index][index_i_center]]){
                            if((*knots[knot_index][index_i_center]).Is_interior()){
                                MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_FOR_JOB, world, &status);
                                work_control[0] = 1;
                                work_control[1] = knot_index;
                                MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
                                //printf("Knot %u is in center [%u %u]\n",knot_index,(*knots[knot_index][0]).index[0],(*knots[knot_index][0]).index[1]);
                                //printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                                (*knots[knot_index][0]).Send_To_Worker(status, world,0);
                                (*knots[knot_index][index_i_center]).Send_To_Worker(status, world,1);
                                center_solved[knots[knot_index][index_i_center]] = true;
                            }
                        }
                }
            }
            for(int worker_index = 0; worker_index < server; worker_index ++){
                Receive_G_B();
            }
            char file_G_i[256],file_G_j[256],file_G[256],file_B_i[256], file_B[256];
            sprintf(file_G_i,"Output/G_i_Det.txt");
            sprintf(file_G_j,"Output/G_j_Det.txt");
            sprintf(file_G,"Output/G_Det.txt");
            sprintf(file_B_i,"Output/B_i_Det.txt");
            sprintf(file_B,"Output/B_Det.txt");
            Write_G_B(file_G_i, file_G_j, file_G, file_B_i, file_B);
            pFile = fopen(debug_fname,"a");
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fprintf(pFile,"Process %d ended its work.\n",myid);
            fclose(pFile);
        } else {
            //Start and end time for the node
            double knot_start, knot_end;
            Eigen::VectorXd X0;
            //G and B storage vectors for each node
            std::vector<double> Bv_temp, G_temp, var_G_temp, bparams_vec(domain_parameters, domain_parameters + 4 * sizeof(domain_parameters[0]));
            std::vector<int>  G_i_temp, G_j_temp, B_i_temp;
            work_control[0] = 1;
            Center supp_center, int_center;
            while(work_control[0] == 1){
                MPI_Send(work_control, 2, MPI_INT, server, ASK_FOR_JOB, world);
                MPI_Recv(work_control, 2, MPI_INT, server, REPLY_WORKER, world, &status);
                if(work_control[0] == 1){
                    knot_start = MPI_Wtime();
                    supp_center.Recieve_From_Server(server,world,0,bvp);
                    int_center.Recieve_From_Server(server,world,1,bvp);
                    Deterministic_Solve(supp_center,int_center, G_i_temp, G_j_temp, G_temp, Bv_temp, B_i_temp);
                        for(unsigned int k = 0;k < G_temp.size(); k ++){
                            G_i.push_back(G_i_temp[k]);
                            G_j.push_back(G_j_temp[k]);
                            G.push_back(G_temp[k]);
                        }
                        for(unsigned int k = 0;k < Bv_temp.size(); k ++){
                            B_i.push_back(B_i_temp[k]);
                            B.push_back(Bv_temp[k]);
                        }
                    Bv_temp.clear();
                    B_i_temp.clear();
                    G_j_temp.clear();
                    G_i_temp.clear();
                    G_temp.clear();
                    knot_end = MPI_Wtime();
                    pFile = fopen("Debug/Node_debug.csv","a");
                    fprintf(pFile,"%d,%f,%f,%f\n",work_control[1],
                    0.0,0.0,(knot_end-knot_start));
                    fclose(pFile);
                }
            }
            pFile = fopen(debug_fname,"a");
            fprintf(pFile,"Process %d is done solving nodes \n", myid);
            Send_G_B();
            fprintf(pFile,"Process %d is done sending G \n", myid);
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fclose(pFile);
        }
        printf("Process %d ended Solve_Interfaces_MC\n",myid);
    }
    void Solve_Interfaces_Semideterministic(BVP bvp, double discretization, unsigned int N_tray){
        double start = MPI_Wtime();
        FILE *pFile;
        if(myid == server){
            pFile = fopen ("Debug/Node_debug.csv","w");
            fprintf(pFile, "index,var_B,APL,time\n");
            fclose(pFile);
            pFile = fopen("Debug/G_var.csv","w");
            if (pFile == NULL) perror("Failed: ");
            fprintf(pFile,"G_i,G_j,var_G_ij\n");
            fclose(pFile);
            for(unsigned int knot_index = 0; knot_index < knots.size(); knot_index ++){
                for(unsigned int index_i_center = 1; index_i_center < knots[knot_index].size(); 
                    index_i_center ++){
                        if(!center_solved[knots[knot_index][index_i_center]]){
                            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_FOR_JOB, world, &status);
                            work_control[0] = 1;
                            work_control[1] = knot_index;
                            MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
                            //printf("Knot %u is in center [%u %u]\n",knot_index,(*knots[knot_index][0]).index[0],(*knots[knot_index][0]).index[1]);
                            //printf("Knot %u is integrated in center [%u %u]\n",knot_index,(*knots[knot_index][1]).index[0],(*knots[knot_index][1]).index[1]);
                            (*knots[knot_index][0]).Send_To_Worker(status, world,0);
                            (*knots[knot_index][index_i_center]).Send_To_Worker(status, world,1);
                            if((*knots[knot_index][index_i_center]).Is_interior())center_solved[knots[knot_index][index_i_center]] = true;
                    }
                }
            }
            for(int worker_index = 0; worker_index < server; worker_index ++){
                Receive_G_B();
            }
            Compute_Solution(bvp);
            pFile = fopen(debug_fname,"a");
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fprintf(pFile,"Process %d ended its work.\n",myid);
            fclose(pFile);
        } else {
            //Start and end time for the node
            double knot_start, knot_end;
            Eigen::VectorXd X0;
            //G and B storage vectors for each node
            double B_temp,B_var_temp;
            std::vector<double> Bv_temp, G_temp, var_G_temp, bparams_vec(domain_parameters, domain_parameters + 4 * sizeof(domain_parameters[0]));
            std::vector<int>  G_i_temp, G_j_temp, B_i_temp;
            //Stencil
            GMSolver solver(bvp ,bparams_vec, discretization, (unsigned int) myid +1);
            work_control[0] = 1;
            Center supp_center, int_center;
            while(work_control[0] == 1){
                MPI_Send(work_control, 2, MPI_INT, server, ASK_FOR_JOB, world);
                MPI_Recv(work_control, 2, MPI_INT, server, REPLY_WORKER, world, &status);
                if(work_control[0] == 1){
                    knot_start = MPI_Wtime();
                    supp_center.Recieve_From_Server(server,world,0,bvp);
                    int_center.Recieve_From_Server(server,world,1,bvp);
                    if(int_center.Is_interior()){
                        Deterministic_Solve(supp_center,int_center, G_i_temp, G_j_temp, G_temp, Bv_temp, B_i_temp);
                        for(unsigned int k = 0;k < G_temp.size(); k ++){
                            G_i.push_back(G_i_temp[k]);
                            G_j.push_back(G_j_temp[k]);
                            G.push_back(G_temp[k]);
                        }
                        for(unsigned int k = 0;k < Bv_temp.size(); k ++){
                            B_i.push_back(B_i_temp[k]);
                            B.push_back(Bv_temp[k]);
                        }
                        Bv_temp.clear();
                        B_i_temp.clear();
                        G_j_temp.clear();
                        G_i_temp.clear();
                        G_temp.clear();
                    }else{
                        G_j_temp.resize(0);
                        G_temp.resize(0);
                        supp_center.Get_knot_position(work_control[1],X0);
                        solver.Solve(supp_center,int_center,work_control[1],N_tray,discretization,
                        G_j_temp,G_temp,var_G_temp,B_temp,B_var_temp);
                        B.push_back(B_temp);
                        //B.push_back(bvp.u.Value(X0,0.0));
                        B_i.push_back((int)work_control[1]);
                        var_B.push_back(B_var_temp);
                        for(unsigned int k = 0;k < G_temp.size(); k ++){
                            G_i.push_back((int)work_control[1]);
                            G_j.push_back(G_j_temp[k]);
                            G.push_back(G_temp[k]);
                        }
                    }
                    knot_end = MPI_Wtime();
                    pFile = fopen("Debug/Node_debug.csv","a");
                    fprintf(pFile,"%d,%f,%f,%f\n",work_control[1],
                    solver.var_B,solver.APL,(knot_end-knot_start));
                    fclose(pFile);
                    pFile = fopen("Debug/G_var.csv","a");
                    for(unsigned int i = 0; i < solver.var_G.size(); i++){
                        fprintf(pFile,"%d,%d,%f\n",work_control[1],
                        G_j_temp[i],solver.var_G[i]);
                    }
                    fclose(pFile);
                }
            }
            pFile = fopen(debug_fname,"a");
            fprintf(pFile,"Process %d is done solving nodes \n", myid);
            Send_G_B();
            fprintf(pFile,"Process %d is done sending G \n", myid);
            double end = MPI_Wtime();
            fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
            fclose(pFile);
        }
        printf("Process %d ended Solve_Interfaces_MC\n",myid);
    }

};
#endif