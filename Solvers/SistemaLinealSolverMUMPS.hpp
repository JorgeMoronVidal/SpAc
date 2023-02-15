#ifndef SISTEMALINEALSOLVER
#define SISTEMALINEALSOLVER
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "../Mesh/SistemaLineal.hpp"
#include <dmumps_c.h>
#define USE_COMM_WORLD -987654
static int _stat = (-1);
void *thread_child(void *arg)
    {
        char buf[256];
        int *pid = (int *)arg;
        unsigned int mem_tmp, mem_min, mem_max;
        double avg_mem;
        avg_mem = 0.0;
        mem_min = (1U << 31) - 1;
        mem_max = 0U;
        int stat0, stat1;
        stat0 = _stat;
        unsigned int count = 0U;
        //  fprintf(stderr, "thread_child forked\n");
        while(_stat != 0) {
            stat1 = _stat;
            if (stat1 == 1) {
                sprintf(buf, "/proc/%d/statm", *pid);
                std::ifstream fin(buf);
                fin >> mem_tmp;
                fin.close();
                if (mem_tmp > mem_max) {
	                mem_max = mem_tmp;
                }
                if (mem_tmp < mem_min) {
	                mem_min = mem_tmp;
                }
                avg_mem += (double)mem_tmp;
                count++;
            }
            if ((stat1 == (-1)) && (stat0 == 1)) {
                fprintf(stderr, 
	            "used memory :min: %14.8e  max: %14.8e avg: %14.8e count: %d\n", 
	            (double)mem_min * 4.0 / (1024.0 * 1024.0),
	            (double)mem_max * 4.0 / (1024.0 * 1024.0),
	            (avg_mem / (double)count) * 4.0 / (1024.0 * 1024.0),
	            count);
                count = 0U;
                avg_mem = 0.0;
                mem_min = (1U << 31) - 1;
                mem_max = 0U;
            }
            stat0 = stat1;
            usleep(1000);
        }
        // fprintf(stderr, "thread_child join\n count = %ld\n", count);
        pthread_exit(arg);
        return (void *)NULL;
    }
class SistemaLinealSolver{
    public:
    std::vector<double> solucion, residuo;
    void SpMV(double *y, double *x, int nrow, int nnz,
	MUMPS_INT *irow, MUMPS_INT *jcol, double *val)
    {
        for (int i = 0; i < nrow; i++) {
                y[i] = 0.0;
        }
        for (int k = 0; k < nnz; k++) {
                y[irow[k] - 1] += val[k] * x[jcol[k] - 1];
        }
    }
    void MUMPS_LU(SistemaLineal & sistema, int decomposer, int num_threads, double eps_pivot, int scaling){
        int n, itmp, jtmp;
        char fname[256], fname1[256];
        int flag;
        std::vector<int> aux_irow, aux_jcol;
        std::vector<double> aux_val;
        MUMPS_INT *irow, *jcol;
        double *val, *coefs;
        int numlevels;
        int minNodes = 64;
        FILE *fp;
        bool upper_flag = true;
        bool kernel_detection_all = false;
        DMUMPS_STRUC_C id;
        int pid = (int)getpid();
        irow = new MUMPS_INT[sistema.nnz];
        jcol = new MUMPS_INT[sistema.nnz];
        val = new double[sistema.nnz];
        sistema.Convierte_COO(aux_irow,aux_jcol,aux_val);
        for(int i = 0; i < sistema.nnz; i++){
            irow[i] = 1 + aux_irow[i] ;
            jcol[i] = 1 + aux_jcol[i] ;
            val[i] = aux_val[i] ;
        }
        getchar();
        aux_irow.resize(0); aux_jcol.resize(0); aux_val.resize(0);
        id.job = (-1); // job init
        //HOST IS INVOLVED INTO  IN THE PARALLEL STEPS OF PARALLELIZATION
        id.par = 1;
        //2 IF USUAL SYMMETRIC, 0 IF NOT
        id.sym = 0;
        id.comm_fortran = USE_COMM_WORLD;
        //INICIALIZATION
        dmumps_c(&id);
        //ANALYSIS
        void* results;
        pthread_attr_t th_attr;
        pthread_t thread;
        pthread_attr_init(&th_attr);
        pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       
        int pthid = pthread_create(&thread, &th_attr, 
			&thread_child,
			(void *)&pid);
        if (pthid != 0) {
            std::cout << "bad thread creation ? " << pid << std::endl;
            exit(0);
        }
        id.job = 1; // symbolic
        id.n = sistema.n_filas;
        id.nz = sistema.nnz;
        id.irn = irow;
        id.jcn = jcol;
        //SIMILAR TO IPARAMS
        id.icntl[0] = 6;
        id.icntl[1] = 6;
        id.icntl[2] = 6;
        id.icntl[3] = 2;
        id.icntl[6] = 3;
        id.icntl[12] = 1;
        id.icntl[15] = num_threads;
        id.icntl[23] = 1;
        dmumps_c(&id);
        id.job = 2; // numeirc
        id.a = val;
        id.icntl[7] = 0;
        id.cntl[2] = (-1.0) * eps_pivot;
        _stat = 1;
        usleep(5000);
        //LU FACTORIZATION
        dmumps_c(&id);
        _stat = (-1);
        usleep(5000);
        double *x = new double[sistema.n_filas];
        double *y = new double[sistema.n_filas];
        double *z = new double[sistema.n_filas];
        for(int i = 0; i < sistema.n_filas; i++){
            y[i] = sistema.B[i];
            z[i] = sistema.B[i];
        }
        id.job = 3; // forward/backward
        id.nrhs = 1;
        id.rhs = y;
        _stat = 1;  
        //SOLUTION IS COMPUTED
        dmumps_c(&id);
        SpMV(x, y, sistema.n_filas, sistema.nnz, irow, jcol, val);
        solucion.resize(sistema.n_filas);
        residuo.resize(sistema.n_filas);
        for(int i = 0; i < sistema.n_filas; i++){
            solucion[i] = y[i];
            residuo[i] = x[i]-z[i];
        }
        _stat = (-1);
        usleep(5000);
        /*delete [] irow;
        delete [] jcol;
        delete [] val;
        delete [] x;
        delete [] y;
        delete [] z;*/
    };
    
};
#endif
