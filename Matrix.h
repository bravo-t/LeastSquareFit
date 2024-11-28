#ifndef _MTX_H_
#define _MTX_H_

#include <vector>

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

Matrix operator*(double scale, const Matrix& m);

#endif
