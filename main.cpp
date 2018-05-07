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
  B.speye(N_grid,N_grid);
  B*= 1-2*alpha;
  B(0,0)=B(N_grid-1,N_grid-1)=1;
  for(int i=N_grid-2; i>0; i--){
    B(i,i+1)=B(i,i-1)=alpha;
  }
}

void plot_vector(arma::Col<double> arma_vec, int N_grid){
  std::vector<double> vec_x(N_grid), vec_y(N_grid);
  for( int i=0; i<N_grid; i++){
    vec_x.at(i)=i;
    vec_y.at(i)=arma_vec(i);
  }
  plt::plot(vec_x,vec_y);
  plt::draw();
  // plt::show();
  plt::pause(0.001);
  // std::cout << "Press enter to continue." << std::endl;
  // getchar();

}

int main(int argc, char *argv[]){
  // Starting timer.
  double start_time=clock();
  double last_time=start_time;

  // Setup of parameters
  double Da0 =4e-7;
  double Db0 =2.0/3.0*Da0;
  double Dc0 =8.0/15.0*Da0;
  double R   =1.0;
  double L   =1.0;
  int T = 1000000;

  // Decide CFL for a, then decide delta_t from this, letting CFL for b and c be smaller.
  double CFL_a = 0.49;
  // Boundary conditions
  double a0  = 1.0;
  double b0  = a0*10.0;
  // Nucleation threshold
  double c0 = 0.03;

  int N_grid = 1000 + 2;
  double delta_x = L/(N_grid-1);
  double delta_t = pow(delta_x,2)*CFL_a/Da0;

  double CFL_b = Db0*delta_t/pow(delta_x,2);
  double CFL_c = Dc0*delta_t/pow(delta_x,2);

  // Setup boundary conditions, aka initial vectors.
  arma::Col<double> a(N_grid);
  a.fill(0);
  a(0)=a0/a0;

  arma::Col<double> b(N_grid);
  b.fill(0);
  b(N_grid-1)=b0/a0;

  // Explicit euler
  // Set up matrices
  arma::SpMat<double> B_a;
  get_euler_explicit_matrix(CFL_a, N_grid, B_a);
  arma::SpMat<double> B_b;
  get_euler_explicit_matrix(CFL_b, N_grid, B_b);

  // Plot
  plot_vector(a,N_grid);
  for(int t=0; t<T; t++){
    a = B_a*a;
    b = B_b*b;
    if(t%10000==0){
      plot_vector(a,N_grid);
      plot_vector(b,N_grid);
    }
  }
  plot_vector(a,N_grid);

  last_time= printMessageTime("Reached end of main.", start_time,last_time);
  return 0;
}
