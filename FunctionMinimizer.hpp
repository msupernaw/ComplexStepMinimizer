/* 
 * File:   FunctionMinimizer.hpp
 * Author: matthewsupernaw
 *
 * Created on September 5, 2014, 10:05 AM
 */

#ifndef FUNCTIONMINIMIZER_HPP
#define	FUNCTIONMINIMIZER_HPP

#include <stdint.h>
#include <complex>
#include <vector>
#include <valarray>
#include <iostream>
#include <iomanip> 


#define H 1e-20


namespace csfm {

    template<class T >
    class FunctionMinimizer {
        typedef std::complex<T> variable;

        std::vector<variable*> parameters;
        std::vector<variable*> active_parameters;
        std::vector<uint32_t> phases;
        T maxc;
        T function_value;
        std::vector<T> gradient;

    public:

        FunctionMinimizer() : maxc(std::numeric_limits<T>::max()) {

        }

        void Register(variable& v, uint32_t phase = 1) {
            parameters.push_back(&v);
            phases.push_back(phase);
        }

        virtual void Initialize() {

        }

        virtual void ObjectiveFunction(variable& f) = 0;

        bool Run() {
            uint32_t max_phase = 0;
            for (int i = 0; i < phases.size(); i++) {
                if (phases[i] > max_phase) {
                    max_phase = phases[i];
                }
            }

            
            bool converged = false;
            for (uint32_t p = 1; p <= max_phase; p++) {
                active_parameters.clear();
                
                for (int i = 0; i < parameters.size(); i++) {
                    if (phases[i] <= p) {
                        active_parameters.push_back(parameters[i]);
                    }
                }
                converged = this->LBFGS(active_parameters, 10000);
            }
            return converged;
            //            std::cout << function_value << "\n";
        }

    private:

        void Print(variable& fv, std::vector<variable*>& p, std::valarray<T>& g) {
            std::cout << "Function Value: " << fv.real() << "\n";
            for (int i = 0; i < p.size(); i++) {
                std::cout << i << " "/* << std::scientific */ << p[i]->real() << " " << g[i] << " | ";
                if (i < p.size() - 1) {
                    std::cout << i + 1 << " " << std::scientific << p[i + 1]->real() << " " << g[i + 1] << " | ";
                    i++;
                    std::cout << std::endl;
                } else {
                    std::cout << std::endl;
                }

            }
            std::cout << std::endl;
        }

        void CallObjectiveFunction(variable& fv) {
            this->ObjectiveFunction(fv);
            //            function_value = fv.real();
        }

        void CallGradient(variable& fv, std::vector<variable*>& p, std::valarray<T>& g) {
            //            g.resize(p.size());
            maxc = std::numeric_limits<T>::min();
            for (int i = 0; i < p.size(); i++) {
                T old = p[i]->imag();
                p[i]->imag(H);
                this->ObjectiveFunction(fv);
                g[i] = fv.imag() / H;
                if (std::fabs(g[i]) > maxc) {
                    maxc = std::fabs(g[i]);
                }
                p[i]->imag(old);
            }
            this->ObjectiveFunction(fv);
        }

        /**
         * Compute the dot product of two vectors.
         * @param a
         * @param b
         * @return 
         */
        const T Dot(const std::valarray<T> &a, const std::valarray<T> &b) {
            T ret = 0;
            for (size_t i = 0; i < a.size(); i++) {
                ret += a[i] * b[i];
            }
            return ret;
        }

        /**
         * Compute the Norm of the vector v.
         *  
         * @param v
         * @return 
         */
        const T Norm(std::valarray<T> &v) {

            T ret = (T) 0.0;
            unsigned int i;
            for (i = 0; i < v.size(); i++) {
                ret += v[i] * v[i];

            }
            return std::sqrt(ret);
        }

        /**
         * returns the a column of a matrix as a std::valarray.
         * @param matrix
         * @param column
         * @return 
         */
        const std::valarray<T> Column(std::valarray<std::valarray<T> > &matrix, size_t column) {

            std::valarray<T> ret(this->active_parameters.size());

            for (int i = 0; i < ret.size(); i++) {
                ret[i] = matrix[i][column];
            }
            return ret;
        }

        /**
         * A modified version of Charles Dubout's LBFGS algorithm released under the 
         * terms of GPL. 
         * http://www.idiap.ch/~cdubout/code/lbfgs.cpp
         * 
         * 
         * @param parameters
         * @param results
         * @param iterations
         * @param tolerance
         * @return 
         */
        bool LBFGS(std::vector<variable* > &parameters, size_t iterations = 10000, T tolerance = (T(1e-4))) {

            size_t max_history = 10000;

            int maxLineSearches_ = 5000;



            std::valarray<T> x(parameters.size());
            std::valarray<T> best(parameters.size());
            //current gradient
            std::valarray<T> g(parameters.size());


            std::valarray<T> ng(x.size());

            //initial evaluation
            variable fx;



            //Call the objective function and collect stats..
            this->CallObjectiveFunction(fx);
            this->function_value = fx.real();
            //            variable nfx(fx);
            //Historical evaluations
            std::valarray<T> px(parameters.size());
            std::valarray<T> pg(parameters.size());
            std::valarray<std::valarray<T> > dxs(std::valarray<T > (max_history), parameters.size());
            std::valarray<std::valarray<T> > dgs(std::valarray<T > (max_history), parameters.size());

            //set parameters
            for (size_t i = 0; i < g.size(); i++) {
                x[i] = parameters[i]->real();
            }

            std::valarray<T> z(parameters.size());


            this->CallGradient(fx, parameters, g);


            T step = 0.1;
            T relative_tolerance;
            T norm_g;
            const size_t nop = parameters.size();
            std::valarray<T> p;
            std::valarray<T>a;
            for (int i = 0; i < iterations; ++i) {


                if ((i % 10) == 0)
                    this->Print(fx, parameters, g);



                if (maxc < tolerance) {
                    this->Print(fx, parameters, g);
                    return true;
                }

                z.resize(g.size());
                for (int i = 0; i < g.size(); i++)
                    z[i] = g[i];


                if (i > 0) {

                    size_t h = std::min<size_t > (i, max_history);
                    size_t end = (i - 1) % h;

                    //update histories
                    for (size_t r = 0; r < nop; r++) {
                        // std::cout<<r<<" "<<g[r] <<" - "<< pg[r]<<"\n";
                        dxs[r][end] = parameters[r]->real() - px[r];
                        dgs[r][end] = g[r] - pg[r];
                    }

                    p.resize(h);
                    a.resize(h);

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end - j + h) % h;
                        p[k] = 1.0 / Dot(Column(dxs, k), Column(dgs, k));

                        a[k] = p[k] * Dot(Column(dxs, k), z);
                        std::valarray<T> tmp = Column(dgs, k);
                        z -= a[k] * Column(dgs, k);
                    }
                    //                    // Scaling of initial Hessian (identity matrix)
                    z *= Dot(Column(dxs, end), Column(dgs, end)) / Dot(Column(dgs, end), Column(dgs, end));

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end + j + 1) % h;
                        const T b = p[k] * Dot(Column(dgs, k), z);
                        z += Column(dxs, k) * (a[k] - b);


                    }

                }//end if(i>0)




                T descent = 0;
                for (size_t j = 0; j < nop; j++) {
                    px[j] = parameters[j]->real();
                    x[j] = px[j];
                    pg[j] = g[j];
                    descent += z[j] * g[j];
                }//end for


                descent *= T(-1.0); // * Dot(z, g);
                if (descent > T(-0.0000000001) * relative_tolerance /* tolerance relative_tolerance*/) {

                    z = g;
                    iterations -= i;
                    i = 0;
                    step = 1.0;
                    descent = -1.0 * Dot(z, g);
                }//end if



                bool down = false;

                int ls;



                best = x;
                for (ls = 0; ls < maxLineSearches_; ++ls) {
                    // Tentative solution, gradient and loss
                    std::valarray<T> nx = x - step * z;

                    for (size_t j = 0; j < nop; j++) {
                        parameters[j]->real(nx[j]);
                        parameters[j]->imag(0);
                    }



                    this->CallObjectiveFunction(fx);

                    if (fx.real() != fx.real()) {

                        for (size_t j = 0; j < nop; j++) {
                            parameters[j]->real(best[j]);
                        }
                        return false;
                    }



                    if (fx.real() <= this->function_value + tolerance * T(0.0001) * step * descent) { // First Wolfe condition

                        this->CallGradient(fx, parameters, ng);
                        if (down || (-1.0 * Dot(z, ng) >= 0.9 * descent)) { // Second Wolfe condition
                            x = nx;
                            g = ng;
                            //                            fx = fx;
                            this->function_value = fx.real();
                            break;
                        } else {
                            step *= 10.0;
                        }
                    } else {
                        step /= 10.0;
                        down = true;
                    }
                }



                if (ls == maxLineSearches_) {
                    for (size_t j = 0; j < nop; j++) {

                        parameters[j]->real(best[j]);
                    }

                    std::cout << "Max line searches!\n";

                    return false;
                }

            }

            return false;
        }
    };



}


#endif	/* FUNCTIONMINIMIZER_HPP */

