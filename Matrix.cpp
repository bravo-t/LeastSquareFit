#include <cassert>
#include <cmath>
#include <cstdio>
#include "Matrix.h"

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
