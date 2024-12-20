#include "LeastSquareFit.h"
#include "BFGS.h"
#include <cmath>
#include <cstdio>

int main()
{
  BFGSMinimizer m;
  BFGSMinimizer::Function bananaFunc = [](const std::vector<double>& x) {
    double x0 = x[0];
    double x1 = x[1];
    return 100 * std::pow(x1 - x0 * x0, 2) + std::pow(1 - x0, 2);
  };
  m.setFunction(bananaFunc);
  m.setInitX({-1.9, 2});
  m.run();
  const std::vector<double>& coor = m.minCoordinate();
  printf("Min coordinate for Rosenbrock function found at (%f, %f)\n", coor[0], coor[1]);
  
  BFGSMinimizer::Function testFunc = [](const std::vector<double>& x) {
    double x0 = x[0];
    double x1 = x[1];
    return x0*x0 - x1*x1;
  };
  BFGSMinimizer m2;
  m2.setFunction(testFunc);
  m2.setInitX({1, 1});
  m2.run();
  const std::vector<double>& coor2 = m2.minCoordinate();
  printf("Min coordinate for x^2-y^2 found at (%f, %f)\n", coor2[0], coor2[1]);
  return 0;
  {
  std::vector<double> x;
  std::vector<double> y;
  /// Data generated with y=1.2*exp(0.5*x)-1.5*x+rand()
  x.push_back(0);  y.push_back(1.273059058);
  x.push_back(0.1);  y.push_back(1.197750036);
  x.push_back(0.2);  y.push_back(1.105699515);
  x.push_back(0.3);  y.push_back(0.95973696);
  x.push_back(0.4);  y.push_back(0.917948706);
  x.push_back(0.5);  y.push_back(0.871964402);
  x.push_back(0.6);  y.push_back(0.747031404);
  x.push_back(0.7);  y.push_back(0.68209238);
  x.push_back(0.8);  y.push_back(0.638827339);
  x.push_back(0.9);  y.push_back(0.582671595);
  x.push_back(1);  y.push_back(0.50457035);
  x.push_back(1.1);  y.push_back(0.456109755);
  x.push_back(1.2);  y.push_back(0.456000863);
  x.push_back(1.3);  y.push_back(0.41537737);
  x.push_back(1.4);  y.push_back(0.33004882);
  x.push_back(1.5);  y.push_back(0.380664613);
  x.push_back(1.6);  y.push_back(0.368607585);
  x.push_back(1.7);  y.push_back(0.326146708);
  x.push_back(1.8);  y.push_back(0.262961273);
  x.push_back(1.9);  y.push_back(0.287638016);
  x.push_back(2);  y.push_back(0.330556729);
  x.push_back(2.1);  y.push_back(0.339977164);
  x.push_back(2.2);  y.push_back(0.367842484);
  x.push_back(2.3);  y.push_back(0.370351028);
  x.push_back(2.4);  y.push_back(0.453984659);
  x.push_back(2.5);  y.push_back(0.456664873);
  x.push_back(2.6);  y.push_back(0.546374062);
  x.push_back(2.7);  y.push_back(0.621052369);
  x.push_back(2.8);  y.push_back(0.757517399);
  x.push_back(2.9);  y.push_back(0.851423913);
  x.push_back(3);  y.push_back(0.900388997);

  LeastSquareFit lsf;
  lsf.addObservationXData(x);
  lsf.addObservationYData(y);
  LeastSquareFit::Function func = [](const std::vector<double>& params, const std::vector<double>& xs) {
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double x = xs[0];
    return a*std::exp(b*x)+c*x;
  };
  lsf.setFunction(func);
  /*
  lsf.addDerivativeFunction([](const std::vector<double>& params, const std::vector<double>& xs) {
    //double a = params[0];
    double b = params[1];
    //double c = params[2];
    double x = xs[0];
    return std::exp(b*x);
  });
  lsf.addDerivativeFunction([](const std::vector<double>& params, const std::vector<double>& xs) {
    double a = params[0];
    double b = params[1];
    //double c = params[2];
    double x = xs[0];
    return a*std::exp(b*x)*x;
  });
  lsf.addDerivativeFunction([](const std::vector<double>& params, const std::vector<double>& xs) {
    //double a = params[0];
    //double b = params[1];
    //double c = params[2];
    double x = xs[0];
    return x;
  });
  */
  lsf.setInitParams({1,1,-1});
  lsf.setFitMethod(LeastSquareFit::Method::BFGS);
  lsf.run();
  std::vector<double> params = lsf.parameters();
  printf("BFGS Fitted data parameters a = %f, b = %f, c = %f, error = %f\n", params[0], params[1], params[2], lsf.standardError());
  lsf.setInitParams({1,1,-1});
  lsf.setFitMethod(LeastSquareFit::Method::GaussNewton);
  lsf.run();
  params = lsf.parameters();
  printf("GN Fitted data parameters a = %f, b = %f, c = %f, error = %f\n", params[0], params[1], params[2], lsf.standardError());
  }
  {
  std::vector<double> x;
  std::vector<double> y;
  x.push_back(1); y.push_back(3.2939);
  x.push_back(2); y.push_back(4.2699);
  x.push_back(4); y.push_back(7.1749);
  x.push_back(5); y.push_back(9.3008);
  x.push_back(8); y.push_back(20.259);
  LeastSquareFit lsf;
  lsf.addObservationXData(x);
  lsf.addObservationYData(y);
  LeastSquareFit::Function func = [](const std::vector<double>& params, const std::vector<double>& xs) {
    double a = params[0];
    double b = params[1];
    double x = xs[0];
    return a*std::exp(b*x);
  };
  lsf.setFunction(func);
  /*lsf.addDerivativeFunction([](const std::vector<double>& params, const std::vector<double>& xs) {
    //double a = params[0];
    double b = params[1];
    //double c = params[2];
    double x = xs[0];
    return std::exp(b*x);
  });
  lsf.addDerivativeFunction([](const std::vector<double>& params, const std::vector<double>& xs) {
    double a = params[0];
    double b = params[1];
    //double c = params[2];
    double x = xs[0];
    return a*std::exp(b*x)*x;
  });*/
  lsf.setInitParams({1, 1});
  lsf.setFitMethod(LeastSquareFit::Method::BFGS);
  lsf.run();
  std::vector<double> params = lsf.parameters();
  printf("BFGS Fitted data parameters a = %f, b = %f, error = %f\n", params[0], params[1], lsf.standardError());
  lsf.setInitParams({1, 1});
  lsf.setFitMethod(LeastSquareFit::Method::GaussNewton);
  lsf.run();
  params = lsf.parameters();
  printf("GN Fitted data parameters a = %f, b = %f, error = %f\n", params[0], params[1], lsf.standardError());
  }
  return 0;
}