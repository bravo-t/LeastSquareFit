#include <cassert>
#include <cmath>
#include <cstdio>
#include "LeastSquareFit.h"

class Matrix {
  public:
    Matrix() = default;
    Matrix(size_t rows, size_t cols)
    : _rows(rows), _cols(cols)
    {
      _data.assign(_rows * _cols, 0);
    }
    Matrix(const std::vector<double>& data)
    : _rows(data.size()), _cols(1), _data(data)
    {}

    ~Matrix() = default;

    size_t width() const { return _cols; }
    size_t height() const { return _rows; }

    void swapRows(size_t row1, size_t row2);
    void PLUDecomp(Matrix& P, Matrix& L, Matrix& U) const;
    Matrix upperTriagSolve(const Matrix& m) const;
    Matrix lowerTriagSolve(const Matrix& m) const;

    void print(const char* info) const;

    double operator()(size_t row, size_t col) const
    {
      return _data[row * _cols + col];
    }
    
    double& operator()(size_t row, size_t col) 
    {
      return _data[row * _cols + col];
    }
    
    double& operator()(size_t row) 
    {
      assert(_cols == 1);
      return _data[row];
    }

    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scale) const;
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator/(const Matrix& other) const;
    Matrix operator/(double div) const;

    double asDouble() const;
    Matrix transpose() const;

    std::vector<double> data() const { return _data; }
    double L2Norm() const;

    static Matrix identity(size_t size)
    {
      Matrix result(size, size);
      for (size_t i=0; i<size; ++i) {
        result(i, i) = 1;
      }
      return result;
    }

  private:
    double rowMultiplySubtract(size_t row1, size_t row2, size_t pivot);

  private:
    size_t _rows = 0;
    size_t _cols = 0;
    std::vector<double> _data;
};

void
Matrix::print(const char* info) const
{
  printf("%s\n", info);
  for (size_t i=0; i<_rows; ++i) {
    for (size_t j=0; j<_cols; ++j) {
      printf("%f ", this->operator()(i, j));
    }
    printf("\n");
  }
}

Matrix
Matrix::upperTriagSolve(const Matrix& m) const
{
  Matrix x(_rows, 1);
  for (int i=_rows-1; i>=0; --i) {
    double sum = 0;
    for (size_t j=i+1; j<_cols; ++j) {
      sum += this->operator()(i, j) * x(j);
    }
    double rem = m(i, 0) - sum;
    x(i, 0) = rem / this->operator()(i, i);
  }
  return x;
}

Matrix
Matrix::lowerTriagSolve(const Matrix& m) const
{
  Matrix x(_rows, 1);
  for (size_t i=0; i<_rows; ++i) {
    double sum = 0;
    for (size_t j=0; j<i; ++j) {
      sum += this->operator()(i, j) * x(j);
    }
    double rem = m(i, 0) - sum;
    x(i, 0) = rem / this->operator()(i, i);
  }
  return x;
}

void
Matrix::swapRows(size_t row1, size_t row2)
{
  for (size_t j=0; j<_cols; ++j) {
    std::swap(_data[row1 * _cols + j], _data[row2 * _cols + j]);
  }
}

double 
Matrix::rowMultiplySubtract(size_t row1, size_t row2, size_t pivot)
{
  double m = this->operator()(row2, pivot) / this->operator()(row1, pivot);
  for (size_t j=0; j<_cols; ++j) {
    this->operator()(row2, j) -= (m * this->operator()(row1, j));
  }
  return m;
}

void 
Matrix::PLUDecomp(Matrix& P, Matrix& L, Matrix& U) const
{
  assert(_rows == _cols);
  L = Matrix::identity(_rows);
  U = *this;
  P = Matrix::identity(_rows);
  for (size_t i=0; i<_rows; ++i) {
    size_t pivot = i;
    double pivotValue = std::abs(U(i, i));
    for (size_t j=i+1; j<_rows; ++j) {
      if (std::abs(U(j, i)) > pivotValue) {
        pivot = j;
        pivotValue = std::abs(U(j, i));
      }
    }
    if (pivot != i) {
      P.swapRows(i, pivot);
      U.swapRows(i, pivot);
    }
  }
  for (size_t i=0; i<_rows; ++i) {
    for (size_t j=i+1; j<_rows; ++j) {
      double m = U.rowMultiplySubtract(i, j, i);
      L(j, i) = m;
    }
  }
}



double
Matrix::L2Norm() const
{
  assert(_rows == 1 || _cols == 1);
  double result = 0;
  for (size_t i=0; i<height(); ++i) {
    for (size_t j=0; j<width(); ++j) {
      result += this->operator()(i, j) * this->operator()(i, j);
    }
  }
  return std::sqrt(result);
}

Matrix 
Matrix::operator*(const Matrix& other) const
{
  assert(_cols == other._rows);
  Matrix result(_rows, other._cols);
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      for (size_t k=0; k<_cols; ++k) {
        result(i, j) += this->operator()(i, k) * other(k, j);
      }
    }
  }
  return result;
}

Matrix 
Matrix::operator*(double scale) const
{
  Matrix result(_rows, _cols);
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      result(i, j) += this->operator()(i, j) * scale;
    }
  }
  return result;
}

Matrix 
operator*(double scale, const Matrix& m)
{
  Matrix result(m.height(), m.width());
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      result(i, j) += m(i, j) * scale;
    }
  }
  return result;
}

Matrix 
Matrix::operator+(const Matrix& other) const
{
  assert(_cols == other._cols && _rows == other._rows);
  Matrix result(_rows, _cols);
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      result(i, j) += this->operator()(i, j) + other(i, j);
    }
  }
  return result;
}

Matrix 
Matrix::operator-(const Matrix& other) const
{
  assert(_cols == other._cols && _rows == other._rows);
  Matrix result(_rows, _cols);
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      result(i, j) += this->operator()(i, j) - other(i, j);
    }
  }
  return result;
}

double
Matrix::asDouble() const 
{
  assert(_cols == 1 && _rows == 1);
  return _data[0];
}

Matrix 
Matrix::operator/(double div) const
{
  Matrix result(_rows, _cols);
  for (size_t i=0; i<result.height(); ++i) {
    for (size_t j=0; j<result.width(); ++j) {
      result(i, j) += this->operator()(i, j) / div;
    }
  }
  return result;
}

Matrix 
Matrix::operator/(const Matrix& other) const
{
  double div = other.asDouble();
  return this->operator/(div);
}

Matrix 
Matrix::transpose() const
{
  Matrix t(_cols, _rows);
  for (size_t i=0; i<_rows; ++i) {
    for (size_t j=0; j<_cols; ++j) {
      t(j, i) = this->operator()(i, j);
    }
  }
  return t;
}

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