
An example of using complex-step differntiation in function minimization. The complex-step derivative approximation is a very convenient way of estimating derivatives numerically. It is a simple and accurate means of finding derivatives of a quantity calculated by an existing algorithm.

This minimizer uses a modified version of Charles Dubout's LBFGS algorithm. 

To use this function minimizer:

```cpp
#include "FunctionMinimizer.hpp"

template<class T>
class MyModel : public csfm::FunctionMinimizer<T> {
    typedef std::complex<T> variable;

    variable x;

public:

    void Initialize() {
        this->Register(x);
    }

    void ObjectiveFunction(variable& f) {
        f = std::sin(x);
    }


};


int main(int argc, char** argv) {
    MyModel<double> fm;
    fm.Initialize();
    if (!fm.Run()) {
        std::cout << "model did not converge!\n";
    } else {

        std::cout << "model converged!\n";
    }

    return 0;
}

```

Output:

Function Value: 0
0 0 1 | 

Function Value: -0.982548
0 -1.3837 0.186011 | 

Function Value: -0.997903
0 -1.50602 0.0647328 | 

Function Value: -0.999745
0 -1.54823 0.0225658 | 

Function Value: -0.999969
0 -1.56293 0.00786798 | 

Function Value: -0.999996
0 -1.56805 0.00274339 | 

Function Value: -1
0 -1.56984 0.000956559 | 

Function Value: -1
0 -1.57046 0.000333532 | 

Function Value: -1
0 -1.57068 0.000116295 | 

Function Value: -1
0 -1.5707 9.41992e-05 | 

model converged!

