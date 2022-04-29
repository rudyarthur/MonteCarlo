#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <functional>


#include "point_generator.hpp"
#include "bins.hpp"
#include "region.hpp"

using std::cout; 
using std::endl;
using std::vector;
using std::function;
using std::string;
using std::min;
using std::max;



class integrate_result{
    public:
        double average;
        double variance;
        double integral;
        double error;
        unsigned calls;

        integrate_result(){
            average = 0;
            integral = 0;
            variance = 0;
            error = 0;
            calls = 0;
        }

        void normalise(double sum, double sum2, double volume, unsigned N){
            average = sum/N;
            integral = average * volume;
            variance = (sum2/N) - (average*average);
            error = sqrt(variance / N) * volume;
            calls = N;
        }

        string result_string(){ 
            std::ostringstream oss;
            //oss << std::setprecision(15) << std::scientific << r.volume * result << " +/- " << r.volume * sqrt(var);
            oss << "I = " << integral << " +/- " << error << " <f> = " << average << " V[f] = " << variance;
            return oss.str();
        }
};

class integrator{
    public:
        point_generator &p; //create random points
        unsigned dim;       //number of dimensions of integral
        integrate_result res;
        double nevals;

        integrator(point_generator &_p) : p(_p) { nevals = 0; }
        virtual integrate_result integrate(function<double(vector<double>&)> f, region &r, unsigned N) = 0;
        string result_string(){ return res.result_string(); }

        integrate_result accumulate(function<double(vector<double>&)> f, region &r, unsigned N){
                double sum = 0;
                double sum2 = 0;
                vector<double> x(r.dim);

                for(unsigned i=0; i<N; ++i){
                    p.generate(x, r);
                    double val = f(x); 
                    ++nevals;
                    sum += val;
                    sum2 += val*val;
                }
                res.normalise(sum, sum2, r.volume, N);
                return res;
        }
};










/*
class integrator{
    public:
        point_generator &p;
        vector<bins> cur_b;
        vector<bins> b;
        unsigned dim, B;
        double result;
        double nevals;
        double subtraction;
        double zerocorrection;

        integrator(point_generator &_p, unsigned _dim, unsigned _B) : p(_p), dim(_dim), B(_B) {
            for(unsigned i=0; i<dim; ++i){ cur_b.push_back( bins(B) ); }
            nevals = 0;
            result = 0;
            subtraction = 0;
            //zerocorrection = 0;
        }
        
        double integrate(function<double(vector<double>&)> f, unsigned N, bool subtract=false){
            vector<double> x(dim); 
            double sum = 0;
            for(unsigned i=0; i<N; ++i){
                p.generate(x);
                double val = f(x); //+zerocorrection;
                val -= (subtract) ? approx(x) : 0;
                ++nevals;
                sum += val;
                for(unsigned j=0; j<dim; ++j){
                    cur_b[j].accumulate(x[j], val);
                }
            }
            result = sum/N;
            return result + subtraction; //- zerocorrection;
        }  



        vector<double> bin_ids(vector<double> &x){
            vector<double> bx(x.size());
            for(unsigned i=0; i<x.size(); ++i){ bx[i] = b[i].bin_id(x[i]); }
            return bx;
        }
        /*
        I = int_0^1 f(r) = int_0^1 g(x) int_0^1 h(y) = Ix Iy

        int_0^1 g(x) = int_b0 g(x) + int_b1 g(x) + int_b2 g(x) + ... 
        = sum_{b in bins} int_b g(x)
        |b| int_b g(x) =  <g>_b 

        The sum in the bin accumulator is
        sx_b = <g>_b int_0^1 h(y) = <g>_b Iy

        r ~ (xb_i, yb_j)
        want <g>_i <h>_j ~ f(r)

        sx_i * sy_j = 
        ( <g>_i Iy ) * (  <h>_j Ix ) = 
        <g>_i <h>_j Ix Iy

        f(r) ~ (sx_i * sy_j)/ ( Ix Iy )
        */
        /*double approx(vector<double> &x){
            vector<double> bx = bin_ids(x);
            double out = 1;
            for(unsigned i=0; i<dim; ++i){ out *= b[i].sum[bx[i]]; }
            return out/pow(result,dim-1);
        }

        /*
        Exact value of approx integral
        (sum_i   sx_i) * (sum_j sy_j)/ ( Ix Iy )
        */
      /*  void initialise_approx(){ 
            b = cur_b;
            for(auto &e : b){ e.normalise(); } 

            /*if(result < 0.1){ //integral is "small"
                zerocorrection = 1; 
                result += zerocorrection; //pretend we are integrating 1+f(x) 
                for(unsigned d=0; d<dim; ++d){
                    for(unsigned i=0; i<b[d].NB; ++i){
                        b[d].sum[i] += zerocorrection;
                    }
                } 
            }*/

           /* subtraction = 1;
            for(unsigned d=0; d<dim; ++d){
                double axis_sum = 0;
                for(unsigned i=0; i<b[d].NB; ++i){axis_sum += b[d].sum[i] * b[d].widths[i];}
                subtraction *= axis_sum;
            }
            subtraction /= pow(result,dim-1);



        }
};*/

#endif