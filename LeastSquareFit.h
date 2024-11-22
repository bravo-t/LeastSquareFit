#ifndef _LSFIT_H_
#define _LSFIT_H_

#include <vector>
#include <functional>

class LeastSquareFit {
  public:
    LeastSquareFit();
    ~LeastSquareFit();

    using FitFunction = std::function<double(const std::vector<double>& params, const std::vector<double>& x)>;
    using FitVecFunction = std::function<std::vector<double>(const std::vector<double>& params, const std::vector<double>& x)>;
    void setFitFunction(const FitFunction& func);
    void setFitFunction(const FitVecFunction& func);
    void addDerivativeFunction(const FitFunction& devFunc);

    void run();
    std::vector<double> parameters() const { return _params; }
    double standardError() const { return _error; }

  private:
    bool useNumericalDerivative() const { return _derivatives.empty(); }
    double numericalDerivative() const;
    double derivative() const;

  private:
    FitVecFunction                   _fitFunction;
    std::vector<FitFunction>         _derivatives;
    std::vector<double>              _params;
    double                           _error;
    std::vector<std::vector<double>> _obsX;
    std::vector<std::vector<double>> _obsY;
};
        
#endif