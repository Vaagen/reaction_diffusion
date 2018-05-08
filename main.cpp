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

void step_func(arma::Col<double>& c, double c0, arma::Col<double>& above_critical_c){
  for(int i=0; i<c.n_elem; i++){
    above_critical_c(i) = c(i) > c0;
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
  // Tube length
  double L   =1.0;
  // Diffusion constants
  double Da0 =4e-7/(L*L);
  double Db0 =2.0/3.0*Da0;
  double Dc0 =8.0/15.0*Da0;
  // Reaction rate a + b -> c
  double R   =1.0/L;
  // Reaction rate for nucleation
  double N1  =R/10.0;
  // Reaction rate for deposition on solid
  double N2  =R/10.0;
  // Simulation time
  int T = 1000000;

  // Decide CFL for a, then decide delta_t from this, letting CFL for b and c be smaller.
  // CFL 0.08 highest possible for explicit scheme for b0=10*a0
  double CFL_a = 0.08;
  // Boundary conditions
  double a0  = 1.0;
  double b0  = a0*10.0;
  // Nucleation threshold
  double c0 = 0.03;
  // Grid size
  int N_grid = 1000 + 2;
  // Parameters dependant of the ones above
  double delta_x = L/(N_grid-1);
  double delta_t = pow(delta_x,2)*CFL_a/Da0;
  int T_steps = T/delta_t;
  // To have same time step b and c must have seperate CFL numbers, we vary CFL_a since it's largest
  double CFL_b = Db0*delta_t/pow(delta_x,2);
  double CFL_c = Dc0*delta_t/pow(delta_x,2);

  // Setup initial vectors.
  arma::Col<double> a(N_grid);
  a.fill(0);
  a(0)=a0/a0;
  arma::Col<double> b(N_grid);
  b.fill(0);
  b(N_grid-1)=b0/a0;
  arma::Col<double> c(N_grid);
  c.fill(0);
  arma::Col<double> s(N_grid);
  s.fill(0);
  // To store reaction terms R*a*b
  arma::Col<double> reaction_ab(N_grid);
  reaction_ab.fill(0);
  // To store reaction terms N1*c*c
  arma::Col<double> reaction_cc(N_grid);
  reaction_cc.fill(0);
  // To check if c > c0
  arma::Col<double> above_critical_c(N_grid);
  above_critical_c.fill(0);
  // To store reaction terms N2*c*s
  arma::Col<double> reaction_cs(N_grid);
  reaction_cs.fill(0);


  // Explicit euler
  // Set up matrices
  arma::SpMat<double> B_a;
  get_euler_explicit_matrix(CFL_a, N_grid, B_a);
  arma::SpMat<double> B_b;
  get_euler_explicit_matrix(CFL_b, N_grid, B_b);
  arma::SpMat<double> B_c;
  get_euler_explicit_matrix(CFL_c, N_grid, B_c);

  std::cout << delta_t << '\n';

  // Propagate time and plot.
  for(int t=0; t<T_steps; t++){
    // a + b -> c
    reaction_ab = delta_t*R*(a%b);
    // Nucleation
    step_func(c, c0, above_critical_c);
    reaction_cc = delta_t*N1*(above_critical_c % (c%c));
    // Deposition on existing solids
    reaction_cs = delta_t*N2*(c%s);

    // Add reaction terms
    a = a - reaction_ab;
    b = b - reaction_ab;
    c = c + reaction_ab - reaction_cc - reaction_cs;
    s = s               + reaction_cc + reaction_cs;

    // Diffusion of gasses
    a = B_a*a;
    b = B_b*b;
    c = B_c*c;

    if(t%10000==0){
      // plt::clf();
      // plot_vector(a,N_grid);
      // plot_vector(b,N_grid);
      plot_vector(s,N_grid);
      // plot_vector(1e6*reaction_ab,N_grid);
      }
  }

  last_time= printMessageTime("Reached end of main.", start_time,last_time);
  return 0;
}
