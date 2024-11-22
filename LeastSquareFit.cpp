#include <cassert>
#include "LeastSquareFit.h"


class Matrix {
  public:
    Matrix(size_t rows, size_t cols)
    : _rows(rows), _cols(cols)
    {
      _data.assign(_rows * _cols, 0);
    }

    ~Matrix();

    size_t width() const { return _cols; }
    size_t height() const { return _rows; }

    double operator()(size_t row, size_t col) const
    {
      return _data[row * _cols + col];
    }
    
    double& operator()(size_t row, size_t col) 
    {
      return _data[row * _cols + col];
    }

    Matrix operator*(const Matrix& other) const;

    std::vector<double> data() const { return _data; }
  private:
    size_t _rows = 0;
    size_t _cols = 0;
    std::vector<double> _data;
};

Matrix Matrix::operator*(const Matrix& other) const
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