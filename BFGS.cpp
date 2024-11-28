#include <cassert>
#include <cmath>
#include <cstdio>
#include "Matrix.h"
#include "BFGS.h"

using Function = BFGSMinimizer::Function;

static double 
numericalDerivative(const Function& func, 
                    const std::vector<double>& x, 
                    size_t i)
{
  double h = 1e-6;
  std::vector<double> x_copy = x;
  x_copy[i] += h;
  double result = (func(x_copy) - func(x)) / h;
  return result;
}

static double 
derivative(const Function& derivativeFunc, 
           const std::vector<double>& x)
{
  return derivativeFunc(x);
}

static Matrix
gradFx(const Function& func, 
       const std::vector<Function>& derivatives,
       const std::vector<double>& x)
{
  bool useNumericalDerivative = false;
  if (derivatives.empty()) {
    useNumericalDerivative = true;
  }
  Matrix result(x.size(), 1);
  for (size_t i=0; i<result.height(); ++i) {
    if (useNumericalDerivative) {
      result(i, 0) = numericalDerivative(func, x, i);
    } else {
      result(i, 0) = derivative(derivatives[i], x);
    }
  }
  return result;
}

void
BFGSMinimizer::setFunction(const Function& func)
{
  _fitFunction = func;
}

void
BFGSMinimizer::addDerivativeFunction(const Function& devFunc)
{
  _derivatives.push_back(devFunc);
}

void 
BFGSMinimizer::setInitX(const std::vector<double>& x)
{
  _x = x;
}

void
BFGSMinimizer::run()
{
  double stepSize = 1;
  Matrix Hk = Matrix::identity(_x.size());
  Matrix xk(_x);
  Matrix gFkp1(_x.size(), 1);
  do {
    Matrix gFk = gradFx(_fitFunction, _derivatives, xk.data());
    Matrix xkp1 = xk - (Hk * gFk * stepSize);
    Matrix sk = xkp1 - xk;
    //xk.print("====xk====");
    //xkp1.print("====xkp1====");
    //sk.print("====sk====");
    gFkp1 = gradFx(_fitFunction, _derivatives, xkp1.data());
    Matrix qk = gFkp1 - gFk;
    Matrix qkT = qk.transpose();
    Matrix rhokInv = qkT * sk;
    double rhok = 1 / rhokInv.asDouble();
    Matrix I = Matrix::identity(_x.size());
    Matrix vk = (I - (rhok * qk * sk.transpose()));
    Hk = (vk.transpose() * Hk * vk) + (rhok * sk * sk.transpose());
    xk = xkp1;
  } while (gFkp1.L2Norm() > 1e-6);
  _x = xk.data();
}
