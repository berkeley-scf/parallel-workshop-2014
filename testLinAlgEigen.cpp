#include <iostream>
#include <Eigen/Dense>
#include <R.h>
#include <Rmath.h>

using namespace Eigen;
using namespace std;


//g++ -o testLinAlgEigen testLinAlgEigen.cpp -I/usr/include/eigen3 -I/usr/share/R/include -lRmath -lR -O3 -Wall -fopenmp  
// Eigen package (libeigen3 on Ubuntu) needs to be installed on your system

int main()
{
  //omp_set_num_threads(8); // ???

  int size = 8000;
  MatrixXd mat(size, size);

  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      mat(i,j) = rnorm(0.0, 1.0);
    }
  }
 
  cout << "done with rnorm" << endl;

  // old form
  /* 
  MatrixXd tmp = mat.transpose() * mat;
  cout << "done crossprod" << endl;

  chol.compute(tmp);
  */

  // form 1
  //  LLT<MatrixXd> chol;
  // chol.compute(mat.transpose() * mat);

  // form 2 from Dirk's JSS paper - I think this allows compiler to better optimize,
  // but this is not using threading (Eigen doc says something about only matrix stuff is threaded so maybe these calcs are not
  LLT<MatrixXd> chol(MatrixXd(size,size).setZero().
                     selfadjointView<Lower>().rankUpdate(mat.adjoint()));
  

  cout << "done chol" << endl;
  
  //  llt.compute(mat.cross(mat));
  //  (mat.cross(mat)).llt() 
}
