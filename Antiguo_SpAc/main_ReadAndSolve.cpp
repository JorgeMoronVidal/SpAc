#include "SpAc.hpp"
#include "poisson_3.hpp"
#include "rectangle.hpp"
int main(int argc, char *argv[]){
    SpAc_solver problem(argc,argv);
    double domain_parameters[4];
    domain_parameters[0] = domain_parameters[1] = -20;
    domain_parameters[2] = domain_parameters[3] = 20;
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, std::string> string_init;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    vector_init["F"] = Equation_F;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping_c);
    bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_theta_RBF);
    //problem.Init(domain_parameters,domain_parameters,12,1.2,41,41,bvp,0.005);
    //problem.Init(domain_parameters,7,1.0,48,bvp,0.005);
    //problem.Solve_Interfaces_Semideterministic(bvp,0.0007,10000);
    char file_G_i_MC[256],file_G_j_MC[256],file_G_MC[256],file_B_i_MC[256], file_B_MC[256],
    file_G_i_Det[256],file_G_j_Det[256],file_G_Det[256],file_B_i_Det[256], file_B_Det[256];
    sprintf(file_G_i_MC,"Output/G_i_MC.txt");
    sprintf(file_G_j_MC,"Output/G_j_MC.txt");
    sprintf(file_G_MC,"Output/G_MC.txt");
    sprintf(file_B_i_MC,"Output/B_i_MC.txt");
    sprintf(file_B_MC,"Output/B_MC.txt");
    sprintf(file_G_i_Det,"Output/G_i_Det.txt");
    sprintf(file_G_j_Det,"Output/G_j_Det.txt");
    sprintf(file_G_Det,"Output/G_Det.txt");
    sprintf(file_B_i_Det,"Output/B_i_Det.txt");
    sprintf(file_B_Det,"Output/B_Det.txt");
    problem.Init(domain_parameters,21,1.2,0.021213,88,bvp,0.005);
    problem.Read_G_B_And_Solve(bvp,file_G_i_MC, file_G_j_MC, file_G_MC, file_B_i_MC, file_B_MC,
    file_G_i_Det, file_G_j_Det, file_G_Det, file_B_i_Det, file_B_Det);
}
