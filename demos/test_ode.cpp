#include <iostream>
#include <fstream> 

#ifdef _WIN32
const char* outpath_crank = "C:\\Users\\lukas\\Documents\\Scicomp2\\myASC-ODE\\build\\output_test_ode_crank.txt";
const char* outpath_implicit = "C:\\Users\\lukas\\Documents\\Scicomp2\\myASC-ODE\\build\\output_test_ode_implicit.txt";
const char* outpath_improved = "C:\\Users\\lukas\\Documents\\Scicomp2\\myASC-ODE\\build\\output_test_ode_improved.txt";
const char* outpath_explicit = "C:\\Users\\lukas\\Documents\\Scicomp2\\myASC-ODE\\build\\output_test_ode_explicit.txt";
#else
const char* outpath_crank = "output_test_ode_crank.txt";
const char* outpath_implicit = "output_test_ode_implicit.txt";
const char* outpath_improved = "output_test_ode_improved.txt";
const char* outpath_explicit = "output_test_ode_explicit.txt";
#endif

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};


int main()
{
  double tend = 40*M_PI;
  int steps = 1000;
  double tau = tend/steps;

  {
    Vector<> y = { 1, 0 };  // initializer list
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
    
    CrankNicolson stepper(rhs);
    // ImplicitEuler stepper(rhs);

    std::ofstream outfile (outpath_crank);
    std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

    for (int i = 0; i < steps; i++)
    {
      stepper.doStep(tau, y);
      std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }

  {
    Vector<> y = { 1, 0 };  // initializer list
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
    
    ImplicitEuler stepper(rhs);

    std::ofstream outfile (outpath_implicit);
    std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;

    for (int i = 0; i < steps; i++)
    {
      stepper.doStep(tau, y);
      std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }

  {
    Vector<> y = { 1, 0 };  // initializer list 
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    ImprovedEuler stepper(rhs);
    std::ofstream outfile (outpath_improved);
    std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    for (int i = 0; i < steps; i++)
    {
      stepper.doStep(tau, y);
      std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }
  
  {
    Vector<> y = { 1, 0 };  // initializer list
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
    ExplicitEuler stepper(rhs);
    std::ofstream outfile (outpath_explicit);
    std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    outfile << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
    for (int i = 0; i < steps; i++)
    {
      stepper.doStep(tau, y);
      std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
      outfile << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
    }
  }

  return 0;
}
