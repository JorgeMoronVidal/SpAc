#include "PseudoespectralCirculoSolver.hpp"
void PseudoespectralCirculoSolver::Inicializa(Eigen::Vector2d centro, double radio, unsigned int N_circunferencia, 
    unsigned int N_radio, double theta_0, BVP Problema)
    {
        assert(N_radio%2 == 1);
        assert(N_circunferencia%2 ==0);
        //Translation of Trefethen cheb.m function
        //Puntos de Chebyshev sobre un diametro del circulo,
        Eigen::VectorXd r;
        //Vector auxiliar para construir la matriz de diferenciación en la direccion radial
        Eigen::VectorXi c;
        //Matriz de diferenciacion en la direccion del diametro del circulo,matriz auxiliar
        Eigen::MatrixXd D,R,dR, C;
        r.resize(N_radio+1);
        c.resize(N_radio+1);
        R.resize((N_radio+1),(N_radio+1));
        D.resize((N_radio+1),(N_radio+1));
        C.resize((N_radio+1),(N_radio+1));
        dR.resize((N_radio+1),(N_radio+1));
        for(unsigned int indice_r = 0; indice_r<N_radio+1; indice_r ++){
            r[indice_r] = radio*cos(M_PI*indice_r/N_radio);
            if(indice_r == 0 || indice_r == N_radio){
                c(indice_r) = 2*pow(-1,indice_r);
            }else{
                c(indice_r) = pow(-1,indice_r);
            }
        }
        for(unsigned int indice_r_1 = 0; indice_r_1<N_radio+1; indice_r_1 ++){
            for(unsigned int indice_r_2 = 0; indice_r_2<N_radio+1; indice_r_2 ++){
                R(indice_r_1,indice_r_2) = r[indice_r_1];
                C(indice_r_1,indice_r_2) = c[indice_r_1]*(1.0/c[indice_r_2]);
            }
        }
        dR = R-R.transpose();
        D = C.cwiseQuotient(dR+Eigen::MatrixXd::Identity((N_radio+1),(N_radio+1)));
        double centinela = 0;
        for(unsigned int fila = 0; fila < D.rows(); fila++){
            centinela =0 ;
            for(unsigned int columna = 0; columna < D.cols(); columna ++){
                centinela += D(fila,columna);
            }
            D(fila,fila) = D(fila,fila)-centinela;
        }
        R.resize(0,0);
        dR.resize(0,0);
        C.resize(0,0);
        c.resize(0);
        unsigned int N2 = (N_radio-1)/2,M2 = N_circunferencia/2;
        Eigen::MatrixXd D2aux,D1,D2,E1,E2;
        D2aux = D*D;
        D1.resize(N2+1,N2+1); D2.resize(N2+1,N2+1);
        E1.resize(N2+1,N2+1); E2.resize(N2+1,N2+1);
        for(unsigned int fila = 0; fila < N2+1; fila++){
            for(unsigned int columna = 0; columna < N2+1; columna ++){
                D1(fila,columna) = D2aux(fila,columna);
                D2(fila,columna) = D2aux(fila, N_radio - columna);
                E1(fila,columna) = D(fila,columna);
                E2(fila,columna) = D(fila, N_radio - columna);
            }
        }
        D2aux.resize(0,0);
        D.resize(0,0);
        double dtheta = 2*M_PI/N_circunferencia;
        Eigen::VectorXd theta, aux_vector;
        Eigen::MatrixXd D2t;
        theta.resize(N_circunferencia);
        aux_vector.resize(N_circunferencia);
        for(unsigned int fila = 0 ; fila < N_circunferencia;fila ++){ 
            theta(fila) = (fila+1)*dtheta + theta_0;
            if(fila == 0){
                aux_vector(fila) = -(M_PI*M_PI/(3*dtheta*dtheta))-(1.0/6);
            }else{
                aux_vector(fila) = 0.5*pow(-1,fila+1)/pow(sin(dtheta*fila*0.5),2);
            }
        }
        D2t.resize(N_circunferencia,N_circunferencia);
        //Matriz de Toepliz
        for(unsigned int fila = 0; fila < N_circunferencia;fila ++){
            for(unsigned int columna = 0; columna < N_circunferencia;columna ++){
                D2t(fila,columna) = aux_vector(columna);
            }
            centinela =  aux_vector(N_circunferencia-1);
            for(unsigned int columna = N_circunferencia-1; columna > 0 ; columna --) aux_vector(columna) = aux_vector(columna-1);
            aux_vector(0) = centinela;
        } 
        R = Eigen::MatrixXd::Zero(N2+1,N2+1);
        //Calculo de las componentes del Laplaciano
        for(unsigned int diagonal = 0; diagonal < N2+1; diagonal ++) R(diagonal,diagonal) = 1.0/r(diagonal);
        Eigen::MatrixXd I_M = Eigen::MatrixXd::Identity(N_circunferencia,N_circunferencia);
        Eigen::MatrixXd L_r1,L_r2,L_t,R2;
        L_r1 = Eigen::MatrixXd::Zero((N2+1)*N_circunferencia,(N2+1)*N_circunferencia);
        L_r2 = Eigen::MatrixXd::Zero((N2+1)*N_circunferencia,(N2+1)*N_circunferencia);
        R2 = R*R;
        //Producto de Kroneker
        for (int i = 0; i < D1.rows(); i++){
            for (int j = 0; j < D1.cols(); j++)
            {
                L_r1.block(i*I_M.rows(), j*I_M.cols(), I_M.rows(), I_M.cols()) = (D1+R*E1)(i,j) * I_M;
            }
        }
        I_M.resize(0,0);
        D1.resize(0,0);
        E1.resize(0,0);
        Eigen::MatrixXd IZZI = Eigen::MatrixXd::Zero(N_circunferencia,N_circunferencia);
        IZZI.block(0,M2,M2,M2) = Eigen::MatrixXd::Identity(M2,M2);
        IZZI.block(M2,0,M2,M2) = Eigen::MatrixXd::Identity(M2,M2);
        //Producto de Kroneker
        for (int i = 0; i < D2.rows(); i++){
            for (int j = 0; j < D2.cols(); j++)
            {
                L_r2.block(i*IZZI.rows(), j*IZZI.cols(), IZZI.rows(), IZZI.cols()) = (D2+R*E2)(i,j) * IZZI;
            }
        }
        IZZI.resize(0,0);
        D2.resize(0,0);
        E2.resize(0,0);
        R.resize(0,0);
        L_t = Eigen::MatrixXd::Zero((N2+1)*N_circunferencia,(N2+1)*N_circunferencia);
        for (int i = 0; i < R2.rows(); i++){
            for (int j = 0; j < R2.cols(); j++)
            {
                L_t.block(i*D2t.rows(), j*D2t.cols(), D2t.rows(), D2t.cols()) = R2(i,j) * D2t;
            }
        }
        R2.resize(0,0);
        D2t.resize(0,0);
        //Matriz del operador Laplaciano
        L_pseudospectral = L_r1 + L_r2 + L_t;
        //Subdomain mesh creation
        interpolador.valores.centro = centro;
        interpolador.valores.r.resize(N_circunferencia,N2+1);
        interpolador.valores.t.resize(N_circunferencia,N2+1);
        interpolador.valores.x.resize(N_circunferencia,N2+1);
        interpolador.valores.y.resize(N_circunferencia,N2+1);
        for(unsigned int fila = 0;fila < N_circunferencia; fila ++){
            for(unsigned int columna = 0; columna < N2+1; columna ++){
                interpolador.valores.r(fila,columna) = r(columna);
                interpolador.valores.t(fila,columna) = theta(fila);
                interpolador.valores.x(fila,columna) = r(columna)*cos(theta(fila));
                interpolador.valores.y(fila,columna) = r(columna)*sin(theta(fila));
            }
        }
        Eigen::VectorXd x_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.x.data(), 
                        interpolador.valores.x.cols()*interpolador.valores.x.rows())),
                        y_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.y.data(), 
                        interpolador.valores.y.cols()*interpolador.valores.y.rows())),
                        r_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.r.data(), 
                        interpolador.valores.r.cols()*interpolador.valores.r.rows())),
                        t_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.t.data(), 
                        interpolador.valores.t.cols()*interpolador.valores.t.rows()));
        Eigen::Vector2d aux;
        for(unsigned int fila = 0; fila < x_v.size(); fila++){
            aux(0) = x_v(fila);
            aux(1) = y_v(fila);
            L_pseudospectral(fila,fila) = L_pseudospectral(fila,fila) + Problema.c.Evalua(aux);
        }
        L_pseudospectral.block(0,0,N_circunferencia,x_v.size()) = Eigen::MatrixXd::Zero(N_circunferencia,x_v.size());
        for(unsigned int fila = 0; fila < N_circunferencia; fila ++){
            L_pseudospectral(fila,fila) = 1.0;
        }
        QR.compute(L_pseudospectral);

}
void PseudoespectralCirculoSolver::Resuelve(BVP Problema){
    Eigen::VectorXd x_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.x.data(), 
                        interpolador.valores.x.cols()*interpolador.valores.x.rows())),
                    y_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.y.data(), 
                        interpolador.valores.y.cols()*interpolador.valores.y.rows())),
                    r_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.r.data(), 
                        interpolador.valores.r.cols()*interpolador.valores.r.rows())),
                    t_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.t.data(), 
                        interpolador.valores.t.cols()*interpolador.valores.t.rows())),
                    solucion;
    rhs_pseudospectral.resize(x_v.size());
    solucion.resize(x_v.size());
    Eigen::Vector2d aux;
    for(unsigned int fila = 0; fila < x_v.size(); fila++){
        aux(0) = x_v(fila);
        aux(1) = y_v(fila);
        if(fila >= interpolador.valores.r.rows()){
            rhs_pseudospectral(fila) = (-1.0)*Problema.f.Evalua(aux);
            //solucion(fila) = Problema.u.Evalua(aux);
        }else{
            rhs_pseudospectral(fila) = Problema.g.Evalua(aux);
            //solucion(fila) = Problema.u.Evalua(aux);
        }
    }
    solucion = QR.solve(rhs_pseudospectral);
    interpolador.valores.z = solucion.reshaped(interpolador.valores.r.rows(),interpolador.valores.r.cols());
}

void PseudoespectralCirculoSolver::Resuelve(BVP Problema, std::vector<Nudo> & nudos, std::vector<double> theta_j){
    Eigen::VectorXd x_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.x.data(), 
                        interpolador.valores.x.cols()*interpolador.valores.x.rows())),
                    y_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.y.data(), 
                        interpolador.valores.y.cols()*interpolador.valores.y.rows())),
                    r_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.r.data(), 
                        interpolador.valores.r.cols()*interpolador.valores.r.rows())),
                    t_v(Eigen::Map<Eigen::VectorXd>(interpolador.valores.t.data(), 
                        interpolador.valores.t.cols()*interpolador.valores.t.rows())),
                    solucion;
    
}