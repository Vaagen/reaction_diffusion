#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
// For matrices and such
#include <armadillo>
// matplotlibcpp to make interactive plots
#include "../matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

// Convenient functions
double printTime(double start_time, double last_time){
  double stop_time =clock();
  double lap_time  = (stop_time-last_time)/double(CLOCKS_PER_SEC);
  double total_time= (stop_time-start_time)/double(CLOCKS_PER_SEC);
  std::cout << "    Time since last update(sec): " << lap_time << std::endl;
  std::cout << "    Total time(sec): " << total_time << std::endl;
  return stop_time;
}

double printMessageTime(std::string message, double start_time, double last_time){
  std::cout << message << std::endl;
  return printTime(start_time,last_time);
}

void get_euler_explicit_matrix(double alpha, int N_grid, arma::SpMat<double>& B){
  B.eye(N_grid,N_grid);
  B*= 1-2*alpha;
  B(0,0)=B(N_grid-1,N_grid-1)=1;
  for(int i=N_grid-2; i>0; i--){
    B(i,i+1)=B(i,i-1)=alpha;
  }
}

int main(int argc, char *argv[]){
  // Starting timer.
  double start_time=clock();
  double last_time=start_time;

  // Setup of parameters
  double Da0 =1;
  double Db0 =1;
  double Dc0 =1;
  double R   =1;
  double L   =1;

  double CLF =1;

  double a0  =1;

  int N_grid = 10 + 2;

  // Setup boundary conditions, aka initial vectors.
  arma::Col<double> a(N_grid);
  a.fill(0);
  a(0)=1;
  a(N_grid-1)=2;

  std::cout << a << '\n';

  // Explicit euler
  // Set up matrices
  double alpha=0.25;
  arma::SpMat<double> B;
  get_euler_explicit_matrix(alpha, N_grid, B);


  // Plot


  last_time= printMessageTime("Reached end of main.", start_time,last_time);
  return 0;
}
