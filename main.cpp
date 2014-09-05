/* 
 * File:   main.cpp
 * Author: matthewsupernaw
 *
 * Created on September 5, 2014, 10:05 AM
 */

#include <cstdlib>
#include <fstream>

#include "FunctionMinimizer.hpp"

using namespace std;

template<class T>
class MyModel : public csfm::FunctionMinimizer<T> {
    typedef std::complex<T> variable;
    std::vector<T> x;
    std::vector<T> y;
    int nobs;
    variable a;
    variable b;
    std::vector<variable> v;

public:

    void Initialize() {
        nobs = 500;
        x.resize(nobs);
        y.resize(nobs);
        
        T aa = 1.90909;
        T bb = 4.0782;
        T init = T(1);
        srand(2014); //ensure the same sequence of random numbers 
        for (int i = 0; i < nobs; i++) {
            v.push_back(variable(T(.1)));
            x[i] = init;
            T secret;
            secret = ((double) rand() / (RAND_MAX)) + 1;
            secret *= 1;
            init += .1;
            y[i] = ((aa * x[i]) / (1 + bb * x[i])) * std::exp(secret);
        }
        a = .00001;
        b = .00001;
        this->Register(a);
        this->Register(b);
        for (int i = 0; i < v.size(); i++) {
            this->Register(v[i], 4);
        }
    }

    void ObjectiveFunction(variable& f) {
        variable norm(0);
        for (int i = 0; i < x.size(); i++) {
            variable tmp = ((a * x[i]) / (T(1) + b * x[i])) * v[i]; 
            norm += (tmp - y[i])*(tmp - y[i]);
        }
        //        
        f = norm; 
    }

    void Dump() {
        std::ofstream report("model_output.csv");
        for (int i = 0; i < x.size(); i++) {
            variable tmp = ((a * x[i]) / (T(1) + b * x[i])) * v[i]; 
            report << x[i] << ", " << y[i] << ", " << tmp.real() << "\n";
        }
    }

};

/*
 * 
 */
int main(int argc, char** argv) {
    MyModel<double> fm;
    fm.Initialize();
    if (!fm.Run()) {
        std::cout << "model did not converge!\n";
    } else {
        fm.Dump();
        std::cout << "model converged!\n";
    }
    
    return 0;
}

