#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>

using namespace std;

extern "C" int dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
extern "C" int dsyrk_(char* uplo, char* trans, int* n, int* k, double* alpha,
                    double* a, int* lda, double* beta, double* c, int* ldc);

// compilation:
// g++ -o testLinAlg testLinAlg.cpp -I/usr/share/R/include -llapack -lblas -lRmath -lR -O3 -Wall
// libRmath.so needs to be installed on the machine you're compiling on; at the moment it is not on all the SCF Linux machines

// Note to CJP: use update-alternatives --config libblas.so and libblas.so.3gf to switch BLAS for testing non-threaded default BLAS

int main(){
  int size = 8000;
  int info = 0;
  char uplo = 'U';
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  double* x = new double[size*size];
  double* C = new double[size*size];
  // vector<double> x[size*size]; // doesn't interface nicely with BLAS
  // vector<double> C[size*size];
  for(int i = 0; i < size*size; i++){
    x[i] = rnorm(0.0, 1.0);
    C[i] = 0.0;
  }
  cout << "done with rnorm" << endl;
  dsyrk_(&uplo, &trans, &size, &size, &alpha, x, &size, &beta, C, &size); 
  cout << "done crossprod" << endl;
  dpotrf_(&uplo,&size,C,&size,&info);
  cout << "done chol" << endl;

  return 0;
}
 

