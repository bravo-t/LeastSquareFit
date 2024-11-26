#ifndef _LSFIT_H_
#define _LSFIT_H_

#include <vector>
#include <functional>

class LeastSquareFit {
  public:
    using Function = std::function<double(const std::vector<double>& params, const std::vector<double>& x)>;
    void setFunction(const Function& func);
    void addDerivativeFunction(const Function& devFunc);

    void setInitParams(const std::vector<double>& params);
    void setObservationData(const std::vector<std::vector<double>>& obsX, 
                            const std::vector<std::vector<double>>& obsY)
    {
      _obsX = obsX;
      _obsY = obsY;
    }

    void addObservationXData(const std::vector<double>& x);
    void addObservationYData(const std::vector<double>& x);

    void run() { GaussNewton(); }
    std::vector<double> parameters() const { return _params; }
    double standardError() const { return _error; }

  private:
    void BFGS();
    void GaussNewton();

  private:
    Function                         _fitFunction;
    std::vector<Function>            _derivatives;
    std::vector<double>              _params;
    double                           _error;
    std::vector<std::vector<double>> _obsX;
    std::vector<std::vector<double>> _obsY;
};
        
#endif