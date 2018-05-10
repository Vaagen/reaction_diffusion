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

void step_func(arma::Col<double>& c, double c0, arma::Col<int>& above_critical_c){
  for(int i=0; i<c.n_elem; i++){
    above_critical_c(i) = c(i) > c0;
  }
}

void update_nucleation_times(int t, int N_grid, arma::Col<int>& t_nucleation,
                              arma::Col<int>& above_critical_c){
  for (int j=0; j<N_grid; ++j){
    if(t_nucleation(j)==0 && above_critical_c(j)>0){
      t_nucleation(j)=t;
    }
  }
}

void diffusion_varying_D(double delta_t, double delta_x,int N_grid, double D0,
          arma::Col<double>& new_a, arma::Col<double> a, arma::Col<double> u){
  // Calculate diffusion term
  u *= D0;
  // Euler scheme with centeral difference
  for(int x=1; x<N_grid-1; ++x){
    new_a(x) = a(x) + 0.25*delta_t/delta_x/delta_x*( u(x+1) - u(x-1) )*( a(x+1) - a(x-1) )
              + u(x)*delta_t/delta_x/delta_x*( a(x+1) -2*a(x) + a(x-1) );
  }
  // Note that boundaries are unchanging
}

void plot_vector(arma::Col<double> arma_vec, int N_grid){
  std::vector<double> vec_x(N_grid), vec_y(N_grid);
  for( int i=0; i<N_grid; i++){
    vec_x.at(i)=i*1.0/N_grid;
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
  R=0.1;
  // Reaction rate for nucleation
  double N1  =R/10.0;
  // Reaction rate for deposition on solid
  double N2  =R/10.0;
  // Simulation time
  int T = 700000;
  // Time step
  // delta_t=0.2 gives CFL_a=0.08, highest possible for explicit scheme with b0=10*a0
  double delta_t = 0.1;
  // Framerate in how many timesteps between each frame
  int framerate = 10000;
  // Boundary conditions
  double a0  = 1.0;
  double b0  = a0*10.0;
  // Nucleation threshold
  double c0 = 0.03;
  // Clogging effect
  double s0 = 0.001;
  // Grid size
  int N_grid = 1002;
  // Parameters dependant of the ones above
  double delta_x = L/(N_grid-1);
  int T_steps = T/delta_t;
  // CFL number must be different for the gasses to get same time step
  double CFL_a = Da0*delta_t/pow(delta_x,2);
  double CFL_b = Db0*delta_t/pow(delta_x,2);
  double CFL_c = Dc0*delta_t/pow(delta_x,2);

  // Inform user of parameters
  std::cout << "Delta t: " << delta_t << '\n';
  std::cout << "# of points: " << N_grid << '\n';

  // Output identifier and path, NB: make folder manually
  std::string PATH="output/run_s13/";
  std::cout << "Output files will be saved to " << PATH << '\n';
  // Save parameters
  std::ofstream paramFile;
  std::string filename = PATH + "parameters.txt";
  paramFile.open (filename);
  paramFile << "delta_t = " << '\t' << delta_t << '\n';
  paramFile << "L = " << '\t' << L << '\n';
  paramFile << "delta_x = " << '\t' << delta_x << '\n';
  paramFile << "N_grid = " << '\t' << N_grid << '\n';
  paramFile << "f_rate = " << '\t' << framerate << '\n';
  paramFile << "a0 = " << '\t' << a0 << '\n';
  paramFile << "b0 = " << '\t' << b0 << '\n';
  paramFile << "c0 = " << '\t' << c0 << '\n';
  paramFile << "s0 = " << '\t' << s0 << '\n';
  paramFile << "Da0 = " << '\t' << Da0 << '\n';
  paramFile << "Db0 = " << '\t' << Db0 << '\n';
  paramFile << "Dc0 = " << '\t' << Dc0 << '\n';
  paramFile << "R = " << '\t' << R << '\n';
  paramFile << "N1 = " << '\t' << N1 << '\n';
  paramFile << "N2 = " << '\t' << N2 << '\n';
  paramFile << "CFL_a = " << '\t' << CFL_a << '\n';
  paramFile << "CFL_b = " << '\t' << CFL_b << '\n';
  paramFile << "CFL_c = " << '\t' << CFL_c << '\n';
  paramFile.close();
  std::cout << "All parameters saved to " << filename <<'\n';

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
  arma::Col<int> above_critical_c(N_grid);
  above_critical_c.fill(0);
  // To store reaction terms N2*c*s
  arma::Col<double> reaction_cs(N_grid);
  reaction_cs.fill(0);
  // To save timesteps when nucleation happens
  arma::Col<int> t_nucleation(N_grid);
  t_nucleation.fill(0);
  // To ease computation of how D changes, u = 1/(1+s/s0)
  arma::Col<double> u;
  u.fill(0);
  // Explicit Euler matrices for diffusion terms
  arma::SpMat<double> B_a;
  get_euler_explicit_matrix(CFL_a, N_grid, B_a);
  arma::SpMat<double> B_b;
  get_euler_explicit_matrix(CFL_b, N_grid, B_b);
  arma::SpMat<double> B_c;
  get_euler_explicit_matrix(CFL_c, N_grid, B_c);

  // Propagate time (T+1 to include last step so that we save final output as well.)
  for(int t=0; t<T_steps+1; t++){
    // a + b -> c, NOTE: % is element wise multiplication in Armadillo
    reaction_ab = delta_t*R*(a%b);
    // For 2a + b -> c
    // reaction_ab = delta_t*R*(a%a%b);
    // Nucleation
    step_func(c, c0, above_critical_c);
    reaction_cc = delta_t*N1*(above_critical_c % (c%c));
    // To include nucleation times.
    update_nucleation_times(t, N_grid, t_nucleation, above_critical_c);
    // Deposition on existing solids
    reaction_cs = delta_t*N2*(c%s);
    // Add reaction terms
    a = a - reaction_ab;
    // a = a - 2*reaction_ab;
    b = b - reaction_ab;
    c = c + reaction_ab - reaction_cc - reaction_cs;
    s = s               + reaction_cc + reaction_cs;
    // Diffusion of gasses
    // Varying D
    u = 1.0/(1.0+s/s0);
    diffusion_varying_D(delta_t, delta_x, N_grid, Da0, a, a, u);
    diffusion_varying_D(delta_t, delta_x, N_grid, Db0, b, b, u);
    diffusion_varying_D(delta_t, delta_x, N_grid, Dc0, c, c, u);
    // Constant D
    // a = B_a*a;
    // b = B_b*b;
    // c = B_c*c;
    // Diffusion and reactions at once. A little less stable.
    // a = B_a*a - reaction_ab;
    // b = B_b*b - reaction_ab;
    // c = B_c*c + reaction_ab - reaction_cc - reaction_cs;

    // Output
    if(t%framerate==0){
      // Live plot
      // plot_vector(a,N_grid);
      // plot_vector(b,N_grid);
      // Output vectors
      filename = PATH + "a_t=" + std::to_string(t) + ".dat";
      a.save(filename,arma::raw_ascii);
      filename = PATH + "b_t=" + std::to_string(t) + ".dat";
      b.save(filename,arma::raw_ascii);
      filename = PATH + "c_t=" + std::to_string(t) + ".dat";
      c.save(filename,arma::raw_ascii);
      filename = PATH + "s_t=" + std::to_string(t) + ".dat";
      s.save(filename,arma::raw_ascii);
      if(t%(framerate*100)==0){
        std::string message = " Passed t-step # " + std::to_string(t);
        last_time = printMessageTime(message, start_time, last_time);
      }
    }
  }

  filename = PATH + "t_nucleation.dat";
  t_nucleation.save(filename,arma::raw_ascii);

  printMessageTime("Reached end of main.", start_time,last_time);
  return 0;
}
