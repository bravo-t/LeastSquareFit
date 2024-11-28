#ifndef _BFGS_H_
#define _BFGS_H_

#include <vector>
#include <functional>

class BFGSMinimizer {
  public:
    using Function = std::function<double(const std::vector<double>& x)>;
    void setFunction(const Function& func);
    void addDerivativeFunction(const Function& devFunc);

    void setInitX(const std::vector<double>& x);

    void run();
    std::vector<double> minCoordinate() const { return _x; }

  private:
    Function                         _fitFunction;
    std::vector<Function>            _derivatives;
    std::vector<double>              _x;
};
        
#endif