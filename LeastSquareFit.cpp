#include <cassert>
#include <cmath>
#include <cstdio>
#include "LeastSquareFit.h"
#include "Matrix.h"


using Function = LeastSquareFit::Function;

void
addObsData(const std::vector<double>& data, 
           std::vector<std::vector<double>>& dst)
{
  if (dst.empty() == false) {
    assert(data.size() == dst.size());
  }
  if (dst.empty()) {
    dst.assign(data.size(), std::vector<double>());
  }
  for (size_t i=0; i<data.size(); ++i) {
    dst[i].push_back(data[i]);
  }
}

void
LeastSquareFit::addObservationYData(const std::vector<double>& data)
{
  _obsY.assign(data.begin(), data.end());
}

void
LeastSquareFit::addObservationXData(const std::vector<double>& data)
{
  addObsData(data, _obsX);
}

static double
residualFx(const Function& func, 
           const std::vector<double>& params,
           const std::vector<std::vector<double>>& obsX, 
           const std::vector<double>& obsY,
           size_t i)
{
  return func(params, obsX[i]) - obsY[i];
}

static Matrix
Fx(const Function& func, 
   const std::vector<double>& params, 
   const std::vector<std::vector<double>>& obsX,
   const std::vector<double>& obsY)
{
  Matrix result(obsX.size(), 1);
  for (size_t i=0; i<obsX.size(); ++i) {
    result(i) = residualFx(func, params, obsX, obsY, i);
  }
  return result;
}

static double 
numericalDerivative(const Function& func, 
                    const std::vector<double>& params, 
                    const std::vector<double>& obsX,
                    size_t i)
{
  double h = 1e-6;
  std::vector<double> params_copy = params;
  params_copy[i] += h;
  double result = (func(params_copy, obsX) - func(params, obsX)) / h;
  return result;
}

static double 
derivative(const Function& derivativeFunc, 
           const std::vector<double>& params, 
           const std::vector<double>& obsX)
{
  return derivativeFunc(params, obsX);
}

static Matrix
delFx(const Function& func, 
      const std::vector<Function>& derivatives,
      const std::vector<double>& params, 
      const std::vector<std::vector<double>>& obsX)
{
  bool useNumericalDerivative = false;
  if (derivatives.empty()) {
    useNumericalDerivative = true;
  }
  Matrix result(params.size(), obsX.size());
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      if (useNumericalDerivative) {
        result(i, j) = numericalDerivative(func, params, obsX[j], i);
      } else {
        result(i, j) = derivative(derivatives[i], params, obsX[j]);
      }
    }
  }
  return result;
}

static Matrix
gradFx(const Function& func, 
       const std::vector<Function>& derivatives,
       const std::vector<double>& params, 
       const std::vector<std::vector<double>>& obsX, 
       const std::vector<double>& obsY)
{
  bool useNumericalDerivative = false;
  if (derivatives.empty()) {
    useNumericalDerivative = true;
  }
  Matrix result(params.size(), 1);
  for (size_t i=0; i<result.height(); ++i) {
    double gradSum = 0;
    for (size_t j=0; j<obsX.size(); ++j) {
      double del;
      if (useNumericalDerivative) {
        del = numericalDerivative(func, params, obsX[j], i);
      } else {
        del = derivative(derivatives[i], params, obsX[j]);
      }
      del *= residualFx(func, params, obsX, obsY, j);
      gradSum += del;
    }
    result(i, 0) = gradSum;
  }
  return result;
}

void
LeastSquareFit::setFunction(const Function& func)
{
  _fitFunction = func;
}

void
LeastSquareFit::addDerivativeFunction(const Function& devFunc)
{
  _derivatives.push_back(devFunc);
}

void 
LeastSquareFit::setInitParams(const std::vector<double>& params)
{
  _params = params;
}

void
LeastSquareFit::BFGS()
{
  double stepSize = 1;
  //Matrix delF = delFx(_fitFunction, _derivatives, _params, _obsX);
  Matrix Hk = Matrix::identity(_params.size());
  Matrix xk(_params);
  Matrix gFkp1(_params.size(), 1);
  do {
    Matrix gFk = gradFx(_fitFunction, _derivatives, xk.data(), _obsX, _obsY);
    Matrix xkp1 = xk - (Hk * gFk * stepSize);
    Matrix sk = xkp1 - xk;
    //xk.print("====xk====");
    //xkp1.print("====xkp1====");
    //sk.print("====sk====");
    gFkp1 = gradFx(_fitFunction, _derivatives, xkp1.data(), _obsX, _obsY);
    Matrix qk = gFkp1 - gFk;
    Matrix qkT = qk.transpose();
    Matrix rhokInv = qkT * sk;
    double rhok = 1 / rhokInv.asDouble();
    Matrix I = Matrix::identity(_params.size());
    Matrix vk = (I - (rhok * qk * sk.transpose()));
    Hk = (vk.transpose() * Hk * vk) + (rhok * sk * sk.transpose());
    xk = xkp1;
  } while (gFkp1.L2Norm() > 1e-6);
  _params = xk.data();
}

Matrix 
solveLinearEquations(const Matrix& A, const Matrix& b)
{
  Matrix P;
  Matrix L;
  Matrix U;
  A.PLUDecomp(P, L, U);
  Matrix PTb = P.transpose() * b;
  Matrix y = L.lowerTriagSolve(PTb);
  Matrix x = U.upperTriagSolve(y);
  return x;
}

void
LeastSquareFit::GaussNewton() 
{
  Matrix param(_params);
  Matrix x;
  do {
    Matrix F = Fx(_fitFunction, param.data(), _obsX, _obsY);
    Matrix delF = delFx(_fitFunction, _derivatives, param.data(), _obsX);
    //F.print("===F===");
    //delF.print("===delF===");
    Matrix A = delF * delF.transpose();
    Matrix b = delF * F * -1;
    //A.print("===A===");
    //b.print("===b===");
    x = solveLinearEquations(A, b);
    //x.print("===x===");
    param = param + x;
    //param.print("===param===");
  } while (x.L2Norm() > 1e-6);
  _params = param.data();
}

double 
calcStdErr(const Function& func, const std::vector<double>& params, 
           const std::vector<std::vector<double>>& obsX,
           const std::vector<double>& obsY)
{
  Matrix F = Fx(func, params, obsX, obsY);
  Matrix e = F.transpose() * F;
  return std::sqrt(e.asDouble() / (obsX.size() - 2));
}

void
LeastSquareFit::run() 
{
  switch (_fitMethod) {
    case Method::BFGS:
      BFGS();
      break;
    case Method::GaussNewton:
      GaussNewton();
      break;
    default:
      break;
  }
  _error = calcStdErr(_fitFunction, _params, _obsX, _obsY);
}