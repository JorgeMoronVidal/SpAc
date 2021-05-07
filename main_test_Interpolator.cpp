#include "SPDDS.hpp"
#include "poisson_3.hpp"
#include "rectangle.hpp"
int main(int argc, char *argv[]){
    SPDDS_solver problem(argc,argv);
    double domain_parameters[4];
    domain_parameters[0] = domain_parameters[1] = -1;
    domain_parameters[2] = domain_parameters[3] = 1;
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, std::string> string_init;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping_c);
    bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_RBF);
    problem.Init(domain_parameters,3,48,bvp,1.0);
    problem.Test_Interpolator_Parallel(bvp);
}