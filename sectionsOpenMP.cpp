#include <iostream>
using namespace std;

// compile with:  g++  -L/usr/local/lib -fopenmp sectionsOpenMP.cpp -o sectionsOpenMP 

int main(){
  
  #pragma omp parallel // starts a new team of threads
  {
    cout << "I'm the 0th chunk.\n";  // should get run by each thread

    #pragma omp barrier // if we include this, all the 0th reports should happen first

    #pragma omp sections // divides the team into sections 
    { 
      // everything herein is run only by a single thread
      #pragma omp section 
      {       cout << "I'm the 1st chunk." << endl; }
      #pragma omp section 
      { 
        cout << "I'm the 2nd chunk." << endl;
        cout << "I'm the 3rd chunk." << endl;
      } 
      #pragma omp section 
      { cout << "I'm the 4th chunk. See ya." << endl; }
    } // implied barrier
  }

  return 0;
}
