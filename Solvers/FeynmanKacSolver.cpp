#include"FeynmanKacSolver.hpp"
void FeynmanKacSolver::Actualiza_Incremento(Eigen::Vector2d & incremento, boost::mt19937 & rng, 
    boost::normal_distribution<double> & normalrng, double sqrth){
        incremento(0) = sqrth*normalrng(rng);
        incremento(1) = sqrth*normalrng(rng);
    }
void FeynmanKacSolver::Paso_Absorbente(Eigen::Vector2d & X, Eigen::Vector2d & N, double & Y, double & Z, double & xi,
    double & t, double ji_t, double h, double sqrth, Eigen::Vector2d & incremento, boost::mt19937 & rng, 
    boost::normal_distribution<double> & normalrng, BVP & bvp , unsigned int & Llamadas_RNG){
        Eigen::Matrix2d sigma_aux = bvp.sigma.Evalua(X);
        Actualiza_Incremento(incremento, rng, normalrng, sqrth);
        Llamadas_RNG += 2;
        Z += bvp.f.Evalua(X)*Y*h + bvp.psi.Evalua(X,N)*Y*ji_t;
        xi += Y *(bvp.F.Evalua(X)).dot(incremento);
        Y += bvp.c.Evalua(X)*Y*h + Y *(bvp.mu.Evalua(X)).dot(incremento) + bvp.varphi.Evalua(X,N) * Y * ji_t;
        X += (bvp.b.Evalua(X)-sigma_aux*bvp.mu.Evalua(X))*h + sigma_aux*incremento;
        t += h;
}
void FeynmanKacSolver::Paso_Reflejante(Eigen::Vector2d & X, Eigen::Vector2d & N, Eigen::Vector2d & Npro, double & Y, double & Z ,
    double & xi, double & t, double & ji_t, double rho, double & d_k,  double h, double  sqrth, Eigen::Vector2d & incremento,
    boost::mt19937 & rng, boost::normal_distribution<double> & normalrng, boost::exponential_distribution<double> & exprng
    ,BVP & bvp ,unsigned int & Llamadas_RNG){
        Eigen::Matrix2d sigma_aux = bvp.sigma.Evalua(X);
        Eigen::Vector2d Xp, Nprop, Np;
        double omega,uc,nu;
        if (d_k > -rho){
            do{
                Actualiza_Incremento(incremento, rng, normalrng, sqrth);
                Xp = X + bvp.b.Evalua(X)*h + sigma_aux*incremento;
                omega =  exprng(rng); //Distribución exponencial con parámetro 1/(2*h)
                Llamadas_RNG += 3;
                uc = (N.transpose()*sigma_aux).dot(incremento) +N.transpose().dot(bvp.b.Evalua(X))*h;
                nu = 0.5 *(uc+sqrt(pow((N.transpose()*sigma_aux).norm(),2.0)*omega+pow(uc,2.0)));
                //d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
                if (d_k < -0.0) d_k = 0.0;
                ji_t = std::max(0.0,nu+d_k);
                Xp = Xp - ji_t*N;
            }while((Xp - Nprop).dot(N)>0.0);
            d_k = bvp.frontera_dominio.Evalua_Distancia(Xp,Nprop,Np);
            if(d_k > 0.0){
                //printf("WARNING: The particle didn't enter in the domain  after Lepingle step.\n");
                Xp = Nprop;
                ji_t += d_k;
                d_k = 0.0;
            }
        }else{
            Actualiza_Incremento(incremento, rng, normalrng, sqrth);
            Llamadas_RNG += 2;
            Xp = X + bvp.b.Evalua(X)*h + sigma_aux*incremento;
            ji_t = 0.0;
            d_k = bvp.frontera_dominio.Evalua_Distancia(Xp,Nprop,Np);
            if(d_k > -0.0){
                Xp = Nprop;
                ji_t = d_k;
                d_k = 0.0;
                //std::cout << ji_t << std::endl;
            }
        }
        Z += bvp.f.Evalua(X)*Y*h + bvp.psi.Evalua(X,N)*Y*ji_t;
        xi +=  Y *bvp.F.Evalua(X).dot(incremento);
        Y += bvp.c.Evalua(X)*Y*h + bvp.varphi.Evalua(X,N) * Y * ji_t;
        X = Xp;
        N = Np;
        Npro = Nprop;
        t += - h;
        ji_t = 0.0;
}
outcome FeynmanKacSolver::Dentro(double & distancia, bool & stoppingbc, bool & Neumannbc, Eigen::Vector2d & X, 
    double &ji_t, double & sqrth, Eigen::Vector2d & Npro, Eigen::Vector2d & N, double Gobet_Constant, BVP & bvp){
    Eigen::Vector2d Npro_aux;
    Eigen::Vector2d N_aux;
    double distancia_aux = bvp.frontera_dominio.Evalua_Distancia(X,Npro_aux,N_aux);
    distancia = bvp.frontera_subdominio.Evalua_Distancia(X,Npro,N);
    stoppingbc = true;
    Eigen::Matrix2d sigma_aux;
    outcome status;
    if(distancia_aux >= distancia){
        distancia = distancia_aux;
        Npro = Npro_aux;
        N = N_aux;
        stoppingbc = bvp.frontera_dominio.Evalua_Absorbente(Npro);
        Neumannbc = bvp.frontera_dominio.Evalua_Neumann(Npro);
        sigma_aux = (distancia < 0.0) ? bvp.sigma.Evalua(X) : bvp.sigma.Evalua(Npro);
        if(stoppingbc){
            status = (distancia < Gobet_Constant*(N.transpose()*sigma_aux).norm()*sqrth) ? dentro : parada_dominio;
        }else{
            status = (distancia <= 0) ? dentro : reflejada;
        }
    }else{
        stoppingbc = true;
        sigma_aux = (distancia < 0.0) ? bvp.sigma.Evalua(X) : bvp.sigma.Evalua(Npro);
        status = (distancia < Gobet_Constant*(N.transpose()*sigma_aux).norm()*sqrth) ? dentro : parada_subdominio;
    }
    switch(status){
        
      case reflejada:
        ji_t = distancia;
        X = Npro;
        break;

      default:
        ji_t = 0.0;
        break;
    }
    return status;
}
FeynmanKacSolver::FeynmanKacSolver(void){
    N = 0;
    for(int i = 0; i < 30; i++) sums[i] = 0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            RNG.resize(omp_get_num_threads());
            normal_dist.resize(omp_get_num_threads());
            exp_dist.resize(omp_get_num_threads());
        }
        int id = omp_get_thread_num();
        RNG[id].discard(id*1E+10);
    }
}
FeynmanKacSolver::FeynmanKacSolver(unsigned int MPIrank){
    N = 0;
    for(int i = 0; i < 30; i++) sums[i] = 0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            RNG.resize(omp_get_num_threads());
            normal_dist.resize(omp_get_num_threads());
            exp_dist.resize(omp_get_num_threads());
        }
        int id = omp_get_thread_num();
        RNG[id].discard(MPIrank*1E+12 + id*1E+10);
    }
};
void FeynmanKacSolver::Simulacion_OMP(Eigen::Vector2d X0, unsigned int numero_trayectorias, double discretizacion_temporal,
    double rho, BVP Problema)
    {   
        N += numero_trayectorias;
        h = discretizacion_temporal;
        sqrth = sqrt(h);
        X_tau.resize(numero_trayectorias);
        Y_tau.resize(numero_trayectorias);
        Z_tau.resize(numero_trayectorias);
        outcome_tau.resize(numero_trayectorias);
        tau.resize(numero_trayectorias);
        xi_tau.resize(numero_trayectorias);
        RNGcallsv.resize(numero_trayectorias);
        threads.resize(numero_trayectorias);
        //Part of the algorithm that is going to happen inside a GPU
        #pragma omp parallel
        {
            int id = omp_get_thread_num();
            Eigen::Vector2d X,normal,proyeccion_normal,incremento;
            double Y,Z,xi,t,ji_t,dist,dist_k;
            bool ccabsorbentes, ccNeumann;
            unsigned int RNGCalls_thread = 0;
            outcome estado_local;
            #pragma omp for
            for(unsigned int n = 0; n < numero_trayectorias; n++){
                X = X0;
                Y = 1;
                Z = 0;
                xi = 0;
                ji_t = 0;
                t = 0;
                RNGCalls_thread = 0;
                Dentro(dist,ccabsorbentes,ccNeumann,X,ji_t,sqrth,proyeccion_normal,normal,-0.5826, Problema);
                threads[n] = id;
                do{
                    if(ccabsorbentes){

                        Paso_Absorbente(X,normal,Y,Z,xi,t,ji_t,h,sqrth, incremento, RNG[id],normal_dist[id],
                        Problema,RNGCalls_thread);
                        RNGCalls_thread += 2;
                    }else{

                            Paso_Reflejante(X,normal,proyeccion_normal,Y,Z,xi,t,ji_t,rho,dist_k,h,sqrth,incremento,
                            RNG[id],normal_dist[id],exp_dist[id],Problema,RNGCalls_thread);
                             
                    }
                    dist_k = dist;
                    estado_local = Dentro(dist,ccabsorbentes,ccNeumann,X,ji_t,sqrth,proyeccion_normal,normal,-0.5826, Problema);
                }while( estado_local != outcome::parada_dominio && estado_local != outcome::parada_subdominio);
                X_tau[n] = proyeccion_normal;
                Y_tau[n] = Y;
                Z_tau[n] = Z;
                tau[n] = t;
                xi_tau[n] = xi;
                RNGcallsv[n] = RNGCalls_thread;
                outcome_tau[n] = estado_local;
                threads[n] = omp_get_thread_num();
            }
            //#pragma omp barrier
        }
    }
void FeynmanKacSolver::Reduce_Analytic(BVP Problema, unsigned int numero_trayectorias){
        #pragma omp parallel
        {
            double score_lineal_nvr_thread;
            double sums_local[30];
            for(int i = 0; i < 30; i++) sums_local[i] = 0.0;
            #pragma omp for
            for(unsigned int n = 0; n<numero_trayectorias; n++){
                score_lineal_nvr_thread = Z_tau[n] + Y_tau[n]*Problema.g.Evalua(X_tau[n]);
                sums_local[ScoreLinear] += score_lineal_nvr_thread;
                sums_local[ScoreLinear2] += score_lineal_nvr_thread* score_lineal_nvr_thread;
                sums_local[XiLinear] += xi_tau[n];
                sums_local[XiLinear2] += xi_tau[n]*xi_tau[n];
                sums_local[ScoreLinearVR] += score_lineal_nvr_thread + xi_tau[n];
                sums_local[ScoreLinearVR2] += pow(score_lineal_nvr_thread + xi_tau[n],2);
                sums_local[XiScoreLinear] += xi_tau[n]*(score_lineal_nvr_thread);
                sums_local[RNGCalls] += RNGcallsv[n];
                sums_local[tauLinear] += tau[n];
            }
            #pragma omp critical
            {
                sums[ScoreLinear] += sums_local[ScoreLinear];
                sums[ScoreLinear2] += sums_local[ScoreLinear2];
                sums[XiLinear] += sums_local[XiLinear];
                sums[XiLinear2] += sums_local[XiLinear2];
                sums[ScoreLinearVR] +=  sums_local[ScoreLinearVR];
                sums[ScoreLinearVR2] += sums_local[ScoreLinearVR2];
                sums[XiScoreLinear] += sums_local[XiScoreLinear];
                sums[RNGCalls] += sums_local[RNGCalls];
                sums[tauLinear] += sums_local[tauLinear];
            }
        }

}

void FeynmanKacSolver::Reduce_Analytic(BVP Problema, unsigned int numero_trayectorias,Nudo & nudo, Interfaz interfaz){
    for(int i = 0; i < interfaz.nudos_circunferencia.size(); i++) nudo.G[interfaz.nudos_circunferencia[i].indice_global] = 0.0;
    nudo.B = 0.0;
    #pragma omp parallel
    {
            double score_lineal_nvr_thread, B_local = 0;
            std::map<int,double> G_local;
            for(int i = 0; i < interfaz.nudos_circunferencia.size(); i++) G_local[interfaz.nudos_circunferencia[i].indice_global] = 0.0;
            double sums_local[30];
            for(int i = 0; i < 30; i++) sums_local[i] = 0.0;
            #pragma omp for
            for(unsigned int n = 0; n<numero_trayectorias; n++){
                switch (outcome_tau[n])
                {
                case outcome::parada_dominio :
                    //std::cout << __FILE__ << __LINE__ << "Parada Dominio" << std::endl;
                    score_lineal_nvr_thread = Z_tau[n] + Y_tau[n]*Problema.g.Evalua(X_tau[n]);
                    B_local += score_lineal_nvr_thread + xi_tau[n];
                    break;
                case outcome::parada_subdominio :
                    //std::cout << __FILE__ << __LINE__ << "Parada Subdominio" << std::endl;
                    score_lineal_nvr_thread = Z_tau[n] + Y_tau[n]*Problema.u.Evalua(X_tau[n]);
                    if(sqrt(pow(X_tau[n](0)-Problema.frontera_subdominio.parametros[1],2)+
                    pow(X_tau[n](1)-Problema.frontera_subdominio.parametros[2],2))<
                    Problema.frontera_subdominio.parametros[0]*0.99) std::cout << __FILE__<<" "<<
                     __LINE__ << " ERROR" << std::endl;
                    interfaz.Update_G(Y_tau[n],X_tau[n],Problema,G_local);
                    B_local += Z_tau[n] + xi_tau[n];
                    break;
                default:
                    std::cout <<__FILE__<<__LINE__<< "ERROR" << std::endl;
                    break;
                }
                sums_local[ScoreLinear] += score_lineal_nvr_thread;
                sums_local[ScoreLinear2] += score_lineal_nvr_thread* score_lineal_nvr_thread;
                sums_local[XiLinear] += xi_tau[n];
                sums_local[XiLinear2] += xi_tau[n]*xi_tau[n];
                sums_local[ScoreLinearVR] += score_lineal_nvr_thread + xi_tau[n];
                sums_local[ScoreLinearVR2] += pow(score_lineal_nvr_thread + xi_tau[n],2);
                sums_local[XiScoreLinear] += xi_tau[n]*(score_lineal_nvr_thread);
                sums_local[RNGCalls] += RNGcallsv[n];
                sums_local[tauLinear] += tau[n];
            }
            #pragma omp critical
            {   
                for(std::map<int,double>::iterator it_G = G_local.begin(); 
                it_G != G_local.end(); it_G++) nudo.G[it_G->first] += it_G->second;
                nudo.B += B_local;
                sums[ScoreLinear] += sums_local[ScoreLinear];
                sums[ScoreLinear2] += sums_local[ScoreLinear2];
                sums[XiLinear] += sums_local[XiLinear];
                sums[XiLinear2] += sums_local[XiLinear2];
                sums[ScoreLinearVR] +=  sums_local[ScoreLinearVR];
                sums[ScoreLinearVR2] += sums_local[ScoreLinearVR2];
                sums[XiScoreLinear] += sums_local[XiScoreLinear];
                sums[RNGCalls] += sums_local[RNGCalls];
                sums[tauLinear] += sums_local[tauLinear];
            }
    }
}

void FeynmanKacSolver::Update(void){
        phi = sums[ScoreLinear]/N;
        phiphi = sums[ScoreLinear2]/N;
        phi_sublinear = sums[ScoreSublinear]/N;
        phiphi_sublinear = sums[ScoreSublinear2]/N;
        phi_VR = sums[ScoreLinearVR]/N;
        phiphi_VR = sums[ScoreLinearVR2]/N;
        xi = sums[XiLinear]/N;
        xixi = sums[XiLinear2]/N;
        xi_sublinear = sums[XiSublinear]/N;
        xixi_sublinear = sums[XiSublinear2]/N;
        xiphi = sums[XiScoreLinear]/N;
        xiphi_sublinear = sums[XiScoreSublinear]/N;
        phi_num = sums[ScoreLinearNum]/N;
        phiphi_num = sums[ScoreLinearNum2]/N;
        phi_sublinear_num = sums[ScoreSublinearNum]/N;
        phiphi_sublinear_num = sums[ScoreSublinearNum2]/N;
        phi_VR_num = sums[ScoreLinearVRNum]/N;
        phiphi_VR_num = sums[ScoreLinearVRNum2]/N;
        phi_sublinearVR_num = sums[ScoreSublinearVRNum]/N;
        phiphi_sublinearVR_num = sums[ScoreSublinearVRNum2]/N;
        xi_num = sums[XiLinearNum]/N;
        xixi_num = sums[XiLinearNum2]/N;
        xi_sublinear_num = sums[XiSublinearNum]/N;
        xixi_sublinear_num = sums[XiSublinearNum2]/N;
        xiphi_num = sums[XiScoreLinearNum]/N;
        xiphi_sublinear_num = sums[XiScoreSublinearNum]/N;
        RNGC = sums[RNGCalls];
        APL = sums[tauLinear]/N;
}
void FeynmanKacSolver::Update(Nudo & nudo){
    for(std::map<int,double>::iterator it_G = nudo.G.begin(); 
                it_G != nudo.G.end(); it_G++){
                    //std::cout << __FILE__ << __LINE__ <<" "<< it_G->first <<  " " << (it_G->second)/N << std::endl;
                    nudo.G[it_G->first] = (it_G->second)/N;
                } 
    //std::cout << __FILE__ << __LINE__ << " " << nudo.B/N << std::endl;
    nudo.B = nudo.B/N;
    //getchar();
    Update();
}
void FeynmanKacSolver::Solve_OMP_Analytic(Eigen::Vector2d X0, unsigned int numero_trayectorias, double discretizacion_temporal,
        double rho, BVP Problema){
        N = 0;
        for(int i = 0; i < 30; i++) sums[i] = 0;
        Simulacion_OMP(X0,numero_trayectorias,discretizacion_temporal,rho,Problema);
        Reduce_Analytic(Problema, numero_trayectorias);
        Update();
}

void FeynmanKacSolver::Solve_OMP_Analytic(Nudo & nudo, unsigned int numero_trayectorias, double discretizacion_temporal,
    double rho, BVP Problema, Interfaz interfaz){
    if(nudo.frontera){
        nudo.B = Problema.g.Evalua(nudo.posicion_cartesiana);
    }else{
        N = 0;
        for(int i = 0; i < 30; i++) sums[i] = 0;
        Simulacion_OMP(nudo.posicion_cartesiana, numero_trayectorias,discretizacion_temporal,rho,Problema);
        Reduce_Analytic(Problema, numero_trayectorias, nudo, interfaz);
        Update(nudo);
    }
}
        