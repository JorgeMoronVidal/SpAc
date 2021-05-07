#include "FKACSolver.hpp"
#include "center.hpp"
class GMSolver: public EMFKAC {
    private:
        /*-dist: Distance to the boundary*/
        double dist, norm;
        /*-Are bc's stopping?*/
        bool stoppingbc;
        bool Inside(){
            bool stoppingbc = bvp.boundary.stop(E_P);
            double dist = bvp.boundary.Dist(params, X, E_P, N);

                if(stoppingbc){

                    if( dist < -0.5826*(N.transpose()*sigma).norm()*sqrth){

                        if (t > 0.0) {

                            status = in;

                        } else {

                            status = time_out;

                        }
                    } else {

                        if (t > 0.0) {

                            status = stop;

                        } else {

                            X = E_P;
                            status = time_out;

                        }

                    }

                }else{
                    if( dist <= -0.5826*(N.transpose()*sigma).norm()*sqrth){

                        if (t > 0.0) {

                            status = in;

                        } else {

                            status = time_out;

                        }
                    } else {
                        if (t > 0.0) {
          
                            status = reflect;

                        } else {

                            X = E_P;
                            status = time_out;

                        }
                    }
                }

            switch(status){

            case stop:
                ji_t = 0.0;
                X = E_P;
                return false;
                break;

            case reflect:
                ji_t = dist;
                X = E_P;
                return true;
                break;

            case time_out:
                ji_t = 0.0;
                t = 0.0;
                return false;
                break;

            default:
                ji_t = 0.0;
                return true;
                break;
            }
        }
        /*Updates statistical quantities*/
    void Update_Stat(double sol_0, double xi){
        sol = sol_0 + xi;
        sums[0] += sol;
        sums[1] += sol*sol;
        sums[2] += sol_0;
        sums[3] += sol_0*sol_0;
        sums[4] += xi;
        sums[5] += xi*xi;
        sums[6] += xi*sol_0;
        sums[7] += pow(sol - sol_a,2);
        N_trayectories ++;
        norm = 1.0/N_trayectories;
        err =fabs(mean-sol_a);
        rerr = err/sol_a;
        mean = sums[0] * norm;
        mse = sums[7] * norm;
        var = sums[1]*norm - mean*mean;
        std = sqrt(norm*var);
    }
    public:
        /*Variance of he G matrix and the B vector*/
        std::vector<double> var_G;
        double var_B, APL;
        /*Class initialization
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process 
        -seed is the RNG seed*/
        GMSolver(void);  
        GMSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed){
                    Configure(boundary_value_problem,boundary_parameters,discretization, seed);
            }

        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        void Solve(Center supp_center, Center int_center, unsigned int knot_index, 
                unsigned int Ntray, double discretization, std::vector<int> & G_j, 
                std::vector<double> & G, std::vector<double> & var_G, double & B,
                double & var_B){
                for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
                Eigen::VectorXd X0;
                if(!supp_center.Get_knot_position(knot_index,X0)) printf("Something went wrong\n");
                supp_center.Print();
                increment.resize(X0.size());
                X.resize(X0.size());
                E_P.resize(X0.size());
                N.resize(X0.size());
                //Boundary of the stencil
                Boundary sten_boundary;
                sten_boundary._init_(Circle, Stopping_c);
                //G and B are emptied
                G_j.clear();
                G.clear();
                var_G.clear();
                double b = 0.0, bb = 0.0;
                B = 0.0;
                APL = 0.0;
                /*Control variable 
                -0 if trayectory still inside
                -1 if trayectory exits by stencil's boundary
                -2 if trayectory exits by problem's boundary*/
                uint16_t control;
                Reset(X0,INFINITY);
                printf("Knot %u : [%f %f]\n",knot_index,X0[0], X0[1]);
                //FILE *exitpoints;
                //exitpoints = fopen("Debug/ExitPoints.csv","w");
                //fprintf(exitpoints,"tray,x,y\n");
                //fclose(exitpoints);
                if(bvp.boundary.Dist(params, X0, E_P, N) >= -1E-06){
                    B = bvp.g.Value(X0,t);
                    G_j.push_back((int)knot_index);
                    G.push_back(0);
                    var_G.push_back(0);
                    var_B = 0;
                    printf("OUT Knot [ %f, %f ] Distance = %f N*sig = %f h = %f",X0[0],X0[1],bvp.boundary.Dist(params, X0, E_P, N),(N.transpose()*bvp.sigma.Value(X0,t)).norm(),h);
                    printf("Params [%f,%f,%f,%f]\n", params[0], params[1], params[2], params[3]);
                } else{
                    double h_cent, dist_boundary, dist_stencil;
                    Eigen::VectorXd normal_boundary, normal_stencil, cp_boundary, cp_stencil;
                    normal_stencil.resize(2); normal_boundary.resize(2);
                    cp_boundary.resize(2);cp_stencil.resize(2);
                    h_cent = discretization;
                    h = std::min(h_cent,
                    pow(sten_boundary.Dist(int_center.parameters, X0, E_P, N)/(
                    (N.transpose()*bvp.sigma.Value(X0,t)).norm()*2*0.5826),2.0));
                    h = std::min(h_cent,
                    pow(bvp.boundary.Dist(params, X0, E_P, N)/(
                    (N.transpose()*bvp.sigma.Value(X0,t)).norm()*2*0.5826),2.0));
                    sqrth = sqrt(h);
                    printf("Knot [ %f, %f ] Distance = %f N*sig = %f h = %f\n",X0[0],X0[1],bvp.boundary.Dist(params, X0, E_P, N),(N.transpose()*bvp.sigma.Value(X0,t)).norm(),h);
                    printf("Ntray %u \n",Ntray);
                    for(unsigned int tray = 0; tray < Ntray; tray ++){
                        Reset(X0,INFINITY);
                        control = 0;
                        do{
                            Increment_Update();
                            N_rngcalls += X0.size();
                            VR_CV_Step();
                            APL += h;
                            if(t <= 0){
                                control = 3;
                            }
                            dist_stencil = sten_boundary.Dist(int_center.parameters, X, cp_stencil, normal_stencil);
                            if( dist_stencil> -0.5826*(normal_stencil.transpose()*sigma).norm()*sqrth){
                                E_P = cp_stencil;
                                control = 1;
                                dist_stencil = sten_boundary.Dist(int_center.parameters, cp_stencil, cp_stencil, normal_stencil);
                                if(fabs(dist_stencil)>1E-05){
                                     printf("dist : %f X: %f %f E_P %f %f N %f %f\n",dist_stencil,X(0),X(1),cp_stencil(0),cp_stencil(1),normal_stencil(0),normal_stencil(1));
                                    getchar();
                                }
                                dist_boundary = bvp.boundary.Dist(params, cp_stencil, cp_boundary, normal_boundary);
                                if(dist_boundary > -0.5826*(normal_boundary.transpose()*sigma).norm()*sqrth){
                                    E_P = cp_boundary;
                                    control = 2;
                                }
                                break;
                            }
                            dist_boundary = bvp.boundary.Dist(params, X, cp_boundary, normal_boundary);
                            if(dist_boundary > -0.5826*(normal_boundary.transpose()*sigma).norm()*sqrth){
                                E_P = cp_boundary;
                                control = 2;
                                break;
                            }
                        } while(control == 0);
                        //exitpoints = fopen("Debug/ExitPoints.csv","a");
                        //fprintf(exitpoints,"%u,%f,%f\n",tray, E_P[0],E_P[1]);
                        //fclose(exitpoints);
                        switch (control) {        
                            case 1:
                                int_center.G_Update(E_P,Y,bvp);
                                b += Z + xi;
                                bb += pow(Z+xi,2.0);
                            break;        
                            case 2:
                                b += Z + xi + Y*bvp.g.Value(E_P,t);
                                bb += pow(Z + xi + Y*bvp.g.Value(E_P,t),2.0);
                            break;
                            case 3:
                                b += Z + xi + Y*bvp.p.Value(X,t);
                                bb += pow(Z + xi + Y*bvp.g.Value(X,t),2.0);
                            break;
                            default : 
                                std::cout << "Something went wrong while solving";
                        }
                    }
                    //system("python3 exitpoint_plot.py");
                    //getchar();
                    int_center.Return_G(Ntray,G,var_G,G_j);
                    for(unsigned int i = 0; i < G.size(); i++) var_G[i] = var_G[i]-G[i]*G[i];
                    B = b/Ntray;
                    var_B = bb/Ntray-pow(B,2.0);
                    APL = APL/(h*Ntray);
                    h = h_cent;
                    sqrth = sqrt(h_cent);
                }
            };
};