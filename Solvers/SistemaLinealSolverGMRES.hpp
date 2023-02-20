#ifndef SISTEMALINEALSOLVERGMRS
#define SISTEMALINEALSOLVERGMRS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <mkl.h>
#include <metis.h>
#include <list>
#include "mpi.h"
#include "../Mesh/SistemaLineal.hpp"
#include <mkl_cblas.h>
#include <mkl_trans.h>
#include "../SpAc/GestorMPI.hpp"
typedef MKL_INT BLAS_INT;
#include <time.h>
typedef struct timespec elapsed_t;
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
typedef struct timespec elapsed_t;
#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec);\
((a).tv_nsec = (b).tv_nsec)

double
blas_dot(const int n, const double* x, const int incx, 
		 const double* y, const int incy)
{
  return cblas_ddot((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
}

void
blas_copy(const int n, const double* x, const int incx, double* y,
	  const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(double));
  }
  else {
    cblas_dcopy((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
  }
}

double blas_l2norm(const int n, double *x, const int incX)
{
  double tmp = blas_dot(n, x, incX, x, incX);
  //std::cout << __FILE__ << " " << __LINE__ << " " << tmp << std::endl;
  return sqrt(tmp); //
}

struct csr_matrix {
  MKL_INT nrow, nnz;
  MKL_INT *ia, *ja;
  double *coefs;
};

void SpMV(const int nrow, const csr_matrix a, std::vector<double> &x, std::vector<double> &y)
{
  const double zero(0.0);
  for (int i = 0; i < nrow; i++) {
    y[i] = zero;
    for (int k = a.ia[i]; k < a.ia[i + 1]; k++) {
      int j = a.ja[k];
      y[i] += a.coefs[k] * x[j];
    }
  }
}

void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
  const double _zero(0.0);
  const double _one(1.0);
    if (dy == _zero) {
        cs = _one;
        sn = _zero;
    } else if (fabs(dy) > fabs(dx)) {
        double temp = dx / dy;
        sn = _one / sqrt( _one + temp*temp );
        cs = temp * sn;
    } else {
        double temp = dy / dx;
        cs = _one / sqrt( _one + temp*temp );
        sn = temp * cs;
    }
}


void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
    double temp  =  cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}

void 
Update(std::vector<double> &x, int k, std::vector<double>* h,
       std::vector<double> &s, std::vector<double> *v)
{
  std::vector<double> y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= h[i][i];
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h[j][i] * y[i];
  }

  for (int j = 0; j <= k; j++) {
    for (int n = 0; n < x.size(); n++) {
      x[n] += v[j][n] * y[j];
    }
  }
}
class SistemaLinealSolverGMRES{
    void get_realtime(elapsed_t *tm)
    {
        clock_gettime(CLOCK_MONOTONIC, tm);
    }

    double convert_time(elapsed_t time1, elapsed_t time0)
    {
    double t;
    t = ((double)time1.tv_sec - 
        (double)time0.tv_sec +
        ((double)time1.tv_nsec - 
	    (double)time0.tv_nsec) / 1.0e+9);
    return t;
    }

    int convert_sec(elapsed_t t)
    {
        return (int)t.tv_sec;
    }

    int convert_microsec(elapsed_t t)
    {
    return (int)(t.tv_nsec / 1.0e+3);
    }
    public:
    std::vector<double> solucion, residuo;
    void Resuelve(SistemaLineal & sistema, GestorMPI & gestor, int max_iter, double tol_eps, int nparts,
    int layer_overlap){
      int itmp, jtmp, ktmp;
      char fname[256], fname1[256];
      char buf[1024];
      MKL_INT nrow = sistema.n_filas, nnz = sistema.nnz, nnz_orig;
      int *ptrows, *indcols;
      MKL_INT *irow, *jcol;
      double *val, *coefs, *tmps;
      double *scalrvec, *scallvec;
      csr_matrix a;
      FILE *fp;
      bool isSym  = false;
      int nexcl = 0;
      int bparts, eparts;
      int mpi_id, mpi_procs;
      bool rlscale = false;
      mpi_id = gestor.id;
      mpi_procs = gestor.numprocesos;
      {
      int jtmp = nparts % mpi_procs;
      int itmp = (nparts - jtmp) / mpi_procs;
      if (mpi_id < jtmp) {
        bparts = (itmp + 1) * mpi_id;
        eparts = bparts + itmp + 1;
      }
      else {
        bparts = (itmp + 1) * jtmp + itmp * (mpi_id - jtmp);
        eparts = bparts + itmp;
      }
      }
      fprintf(stderr, "%s %d : %d %d %d\n", __FILE__, __LINE__,
	    mpi_id, bparts, eparts);
      irow = new MKL_INT[nnz];
      jcol = new MKL_INT[nnz];
      val = new double[nnz];
      clock_t t0_cpu, t1_cpu, t2_cpu;
      elapsed_t t0_elapsed, t1_elapsed, t2_elapsed;
      elapsed_t t3_elapsed, t4_elapsed;
      double mpi_elapsed;
      mpi_elapsed = 0.0;
      t0_cpu = clock();
      get_realtime(&t0_elapsed);
      a.nrow = nrow;
      a.nnz = nnz;
      a.ia = new MKL_INT[nrow + 1];
      a.ja = new MKL_INT[nnz];
      a.coefs = new double[nnz];
      for(int i = 0; i < nrow +1; i++){
        a.ia[i] = sistema.G_i[i];
      }
      for(int i = 0; i < nnz; i++){
        a.ja[i] = sistema.G_j[i];
        a.coefs[i] = sistema.G_ij[i];
      }
      idx_t *xadj, *adjcy, *part;
      xadj = new idx_t[nrow + 1];
      adjcy = new idx_t[nnz - nrow];
      part = new idx_t[nrow];
      int *itmps0, *itmps1;
      itmps0 = new int[nrow];
      itmps1 = new int[nrow];
      //std::cout <<" Gestor id " <<  gestor.id << " nnz " << nnz << " nrows " << nrow <<  std::endl;
      for (int i = 0; i < nrow; i++) {
        itmps0[i] = 0;
        for (int k = a.ia[i]; k < a.ia[i + 1]; k++) {
          if (a.ja[k] != i) {
	          itmps0[i]++;
          }
        }
        //printf("i %d  imps0[i] %d\n", i, itmps0[i]);
      }
      xadj[0] = 0;
      for (int i = 0; i < nrow; i++) {
        xadj[i + 1] = xadj[i] + itmps0[i];
        //printf("i %d  xadj[i] %d\n", i, xadj[i]);
      }
      for (int i = 0; i < nrow; i++) {
        itmps0[i] = xadj[i];
        //printf("i %d  imps0[i] %d\n", i, itmps0[i]);
      }
     
      for (int i = 0; i < nrow; i++) {
        for (int k = a.ia[i]; k < a.ia[i + 1]; k++) {
          if (a.ja[k] != i) {
            //printf("nnz %d  i %d  nrow %d  imps0[i] %d k %d \n",nnz,i,nrow, itmps0[i], k);
	          adjcy[itmps0[i]] = a.ja[k];
	          itmps0[i]++;
          }
        }
      }
      int ncon = 1;
      idx_t options[METIS_NOPTIONS] = {0} ;
      int objval;
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_NUMBERING] = 0;
      options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
      fprintf(stderr, "nparts = %d\n", nparts);
      METIS_PartGraphRecursive(&nrow, &ncon, xadj, adjcy, NULL, NULL, NULL, &nparts,
			  NULL, NULL, options, &objval, part);
      fprintf(stderr, "\n");
      delete [] xadj;
      delete [] adjcy;
      int **mask, *mask_work;
      mask_work = new int[nrow * nparts];
      mask = new int*[nrow];
      for (int i = 0; i < nparts; i++) {
        mask[i] = &mask_work[i * nrow];
      }
      memset(mask_work, 0, sizeof(int) * nrow * nparts); // zero clear
      // 1st layer overlap
      for (int i = 0; i < nrow; i++) {
        mask[part[i]][i] = 1;    
      }
      for (int ll = 0; ll < layer_overlap; ll++) {
        for (int n = 0; n < nparts; n++) {
          for (int i = 0; i < nrow; i++) {
	          itmps0[i] = mask[n][i];
          }
          for (int i = 0; i < nrow; i++) {
	          if (itmps0[i] == 1) {
	            for (int k = a.ia[i]; k < a.ia[i + 1]; k++) {
	              mask[n][a.ja[k]] = 1;
	            }
	          }
          }
        }
      }
      csr_matrix *aa;
      //solvers.resize(nparts);
      aa = new csr_matrix[nparts];
      std::vector<int> *submatrix_indx;
      submatrix_indx = new std::vector<int>[nparts];
      double *partition_unity;
      partition_unity = new double[nrow];
      for (int i = 0; i < nrow; i++) {
        partition_unity[i] = 0.0;
      }
      for (int n = 0; n < nparts; n++) {
        for (int i = 0; i < nrow; i++) {
          if (mask[n][i] == 1) {
	          submatrix_indx[n].push_back(i);
	          partition_unity[i] += 1.0;
          }
        }
      }
      for (int i = 0; i < nrow; i++) {
        partition_unity[i] = 1.0 / partition_unity[i];
      }
      for (int n = 0; n < nparts; n++) {
      aa[n].nrow =submatrix_indx[n].size();
      aa[n].ia = new int[aa[n].nrow + 1];

      // inverse index to generate sub sparse matrix
      for (int i = 0; i < aa[n].nrow; i++) {
        itmps1[submatrix_indx[n][i]] = i;
      }
      for (int i = 0; i < nrow; i++) {
        itmps0[i] = 0;
      }
      for (int i = 0; i < aa[n].nrow; i++) {
        int ii = submatrix_indx[n][i];
        for (int k = a.ia[ii]; k < a.ia[ii + 1]; k++) {
	          if (mask[n][a.ja[k]] == 1) {
	            itmps0[i]++;
	          }
        }
      }
      aa[n].ia[0] = 0;
      for (int i = 0; i < aa[n].nrow; i++) {
        aa[n].ia[i + 1] = aa[n].ia[i] + itmps0[i];
      }
      for (int i = 0; i < aa[n].nrow; i++) {
        itmps0[i] = aa[n].ia[i];
      }
      aa[n].nnz = aa[n].ia[aa[n].nrow];
      aa[n].ja = new int[aa[n].nnz];
      aa[n].coefs = new double[aa[n].nnz];
      for (int i = 0; i < aa[n].nrow; i++) {
          int ii = submatrix_indx[n][i];
          for (int k = a.ia[ii]; k < a.ia[ii + 1]; k++) {
	          if (isSym) {
	            if (mask[n][a.ja[k]] == 1 && a.ja[k] >= ii) {
	              aa[n].ja[itmps0[i]] = itmps1[a.ja[k]];
	              aa[n].coefs[itmps0[i]] = a.coefs[k];
	              itmps0[i]++;
	            }
	          } else {
	            if (mask[n][a.ja[k]] == 1) {
	              aa[n].ja[itmps0[i]] = itmps1[a.ja[k]];
	              aa[n].coefs[itmps0[i]] = a.coefs[k];
	              itmps0[i]++;
	            }
	          }
          }
        }
      }
      delete [] mask;
      delete [] mask_work;
      char buff[256];
      for (int n = bparts; n < eparts; n++) {
        sprintf(buff,"Output/Debug/csr_%d.txt",n);
        std::ofstream ff(buff); 
        for (int i = 0; i < aa[n].nrow; i++) {
  	      for(int k = aa[n].ia[i]; k < aa[n].ia[i+1]; k ++){
  		      ff << i+1 << "\t" << aa[n].ja[k]+1 << "\t" << aa[n].coefs[k] << std::endl; 
  	      }
        }
        ff.close();
      }
      DMUMPS_STRUC_C *id = new DMUMPS_STRUC_C[nparts];
      /* When compiling with -DINTSIZE64, MUMPS_INT is 64-bit but MPI
      ilp64 versions may still require standard int for C interface. */
      /* MUMPS_INT myid, ierr; */
      MPI_Comm my_mpi_comm;
      int mpicolor = mpi_id;
      int mpikey = 0;
      MPI_Comm_split(MPI_COMM_WORLD, mpicolor, mpikey, &my_mpi_comm);
      fprintf(stderr, "%s %d :%d %d\n", __FILE__, __LINE__, mpi_id, my_mpi_comm);
      for (int n = bparts; n < eparts; n++) {
        id[n].comm_fortran= (MUMPS_INT) MPI_Comm_c2f(my_mpi_comm);
        id[n].par = 1;
        id[n].sym = 0;
        id[n].job = JOB_INIT;
        dmumps_c(&id[n]);
      }
      double **x, **b;
      x = new double*[nparts];
      b = new double*[nparts];
      
      for (int n = 0; n < nparts; n++) {
        x[n] = new double[aa[n].nrow];
        b[n] = new double[aa[n].nrow];
      }
      for (int n = bparts; n < eparts; n++) {
        std::cout << __FILE__ <<" " << __LINE__ << " " << mpi_id <<" "<< n<< std::endl;
        id[n].icntl[0] = 6;
        id[n].icntl[1] = 6;
        id[n].icntl[2] = 6;
        id[n].icntl[3] = 2;
        id[n].icntl[6] = 3; // compute a symmetric permutation
        id[n].icntl[12] = 1; // controls parallelism of root node : no scalapack
        id[n].icntl[23] = 0; // controls detection of "null pivot rows"
        id[n].n   = aa[n].nrow;
        id[n].nnz = aa[n].nnz;
        id[n].irn = new MUMPS_INT[aa[n].nnz];
        id[n].jcn = new MUMPS_INT[aa[n].nnz];      
        id[n].a = new double[aa[n].nnz];
        id[n].nrhs = 1;
        id[n].rhs = new double[aa[n].nrow];
        int kk = 0; 
        for (int i = 0; i < aa[n].nrow; i++) {
          for (int k = aa[n].ia[i]; k < aa[n].ia[i + 1]; k++) {
	          const int j = aa[n].ja[k];
	          id[n].irn[kk] = i + 1;
	          id[n].jcn[kk] = j + 1;
	          id[n].a[kk] = aa[n].coefs[k];
	          kk++;
          }
        }
        /* Call the MUMPS package (analyse, factorization and solve). */
        id[n].job = 1;
        dmumps_c(&id[n]);
        id[n].job = 2;
        dmumps_c(&id[n]);
      } // loop : n
      //dfgmres_init (&nrow, sol, rhs, &RCI_request, ipar, dpar, tmp);
      std::vector<double> s, cs(max_iter + 1), sn(max_iter + 1), w(nrow);
      std::vector<double> * v = new std::vector<double>[max_iter + 1];
      std::vector<double> sol(nrow), rhs(nrow), exact(nrow);
  
      for (int i = 0; i < max_iter + 1; i++) {
        v[i].resize(nrow);
      }
      std::vector<double> *Hessenberg;
      Hessenberg = new std::vector<double>[max_iter + 1];
      for (int i = 0; i < max_iter + 1; i++) {
        Hessenberg[i].resize(max_iter);
      }
      tmps = new double[nrow];
      for(int i = 0; i < sistema.n_filas; i ++){
          rhs[i] = sistema.B[i];
      }
      for (int i = 0; i < nrow; i++) {
        tmps[i] = 0.0;
      }
      for (int n = bparts; n < eparts; n++) {
        id[n].icntl[0] = -1;
        id[n].icntl[1] = -1;
        id[n].icntl[2] = -1;
        id[n].job = 3;
        for (int i = 0; i < submatrix_indx[n].size(); i++) {
          id[n].rhs[i] = rhs[submatrix_indx[n][i]];
        }
        dmumps_c(&id[n]);
        for (int i = 0; i < submatrix_indx[n].size(); i++) {
        int ii = submatrix_indx[n][i];
          tmps[ii] += id[n].rhs[i] * partition_unity[ii];
        }
      }
      get_realtime(&t3_elapsed);
      MPI_Allreduce(tmps, &rhs[0], nrow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      get_realtime(&t4_elapsed);
      mpi_elapsed += convert_time(t4_elapsed, t3_elapsed);
      for (int i = 0; i < nrow; i++) {
        sol[i] = 0.0;
      }
      t1_cpu = clock();
      get_realtime(&t1_elapsed);
      double beta = blas_l2norm(nrow, &rhs[0], 1);
      std::vector<double> r(nrow);
      blas_copy(nrow, &rhs[0], 1, &r[0], 1);
      for (int n = 0; n < nrow; n++) {
        v[0][n] = r[n] / beta;
      }
      s.resize(max_iter + 1, 0.0);
      s[0] = beta;
      if (mpi_id == gestor.servidor) {
        fprintf(stderr, "beta = %.12e\n", beta);
      }
      int m;
      bool flag_conv = false;
      //GMRES LOOP
      for (m = 0; m < max_iter; m++){
        SpMV(nrow, a, v[m], w);
      for (int i = 0; i < nrow; i++) {
	      tmps[i] = 0.0;
      }
      //Load Balance
      for (int n = bparts; n < eparts; n++) {
        id[n].job = 3;
        id[n].icntl[0] = -1;
        id[n].icntl[1] = -1;
        id[n].icntl[2] = -1;
        for (int i = 0; i < submatrix_indx[n].size(); i++) {
	        id[n].rhs[i] = w[submatrix_indx[n][i]];
        }
        dmumps_c(&id[n]);
        for (int i = 0; i < submatrix_indx[n].size(); i++) {
	      int ii = submatrix_indx[n][i];
	      tmps[ii] += id[n].rhs[i] * partition_unity[ii];
        }
      }
      get_realtime(&t3_elapsed);
      MPI_Allreduce(tmps, &w[0], nrow, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      get_realtime(&t4_elapsed);
      mpi_elapsed += convert_time(t4_elapsed, t3_elapsed);
      //      w[i] = tmps[i];    //    for (int i = 0; i < nrow; i++) {
      for (int k = 0; k <= m; k++) {
        Hessenberg[k][m] = blas_dot(nrow, &w[0], 1, &v[k][0], 1);
        for (int i = 0; i < nrow; i++) {
	        w[i] -= Hessenberg[k][m] * v[k][i];
        }
      }
      double htmp = blas_l2norm(nrow, &w[0], 1);
      Hessenberg[m + 1][m] = htmp;
      //    fprintf(stderr, "%d %g\n", m, htmp);
      htmp = 1.0 / htmp;
      for (int i = 0; i < nrow; i++) {
        v[m + 1][i] = w[i] * htmp;
      }
      for (int k = 0; k < m; k++) {
        ApplyPlaneRotation(Hessenberg[k][m], Hessenberg[k + 1][m],
			  cs[k], sn[k]);
      }
    
      GeneratePlaneRotation(Hessenberg[m][m], Hessenberg[m +1][m],
				    cs[m], sn[m]);
      ApplyPlaneRotation(Hessenberg[m][m], Hessenberg[m + 1][m],
			      cs[m], sn[m]);
      ApplyPlaneRotation(s[m], s[m + 1], cs[m], sn[m]);
      if (mpi_id == gestor.servidor) {
        fprintf(stderr, "%d %.12e\n", m, s[m + 1]);
      }
      if ( (fabs(s[m + 1]) / beta ) < tol_eps) {
        flag_conv = true;
        break;
      }
    } // loop : m
    if (mpi_id == gestor.servidor) {
      fprintf(stderr, "\n");
    }
    Update(sol, flag_conv ? m : (m - 1), Hessenberg, s, v);
    t2_cpu = clock();
    get_realtime(&t2_elapsed);

    if (mpi_id == gestor.servidor) {  
      fprintf(stderr, "%s %d %.4e %.4e\n",  __FILE__, __LINE__,
	      (double)(t2_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t2_elapsed, t0_elapsed));
    
      fprintf(stderr, "%s %d %.4e %.4e\n",  __FILE__, __LINE__,
	      (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t1_elapsed, t0_elapsed));
    
      fprintf(stderr, "%s %d %.4e %.4e\n",  __FILE__, __LINE__,
	      (double)(t2_cpu - t1_cpu) / (double)CLOCKS_PER_SEC,
	      convert_time(t2_elapsed, t1_elapsed));
      fprintf(stderr, "%s %d MPI %.4e\n",  __FILE__, __LINE__,
	      mpi_elapsed);
      std::ofstream file_sol("Output/solution_GMRES.txt");
      file_sol.precision(12);
      file_sol.setf(std::ios::fixed, std::ios::scientific);
      for (int i = 0; i < nrow; i++) {
        file_sol << sol[i] << std::endl;
        sistema.u[i] = sol[i];
      }
      file_sol.close();
  }
  }
};
#endif