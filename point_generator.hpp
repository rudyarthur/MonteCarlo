#ifndef POINT_GENERATOR_H
#define POINT_GENERATOR_H

#include <vector>
#include "sobol/sobol.hpp"
#include "pcg-cpp-0.98/include/pcg_random.hpp"
#include <random>
#include <string>
#include <stdlib.h>
using std::vector;
using std::string;

#include "region.hpp"

class point_generator {
    public:
        string method_name;
        long long int seed;
        std::uniform_real_distribution<double> unif;

        point_generator(long long int _seed) : unif(0,1){
            seed = _seed;
        }
        virtual double random() = 0;
        virtual void generate(vector<double> &x, region &r) = 0;
};

class pcg_generator : public point_generator {
    public:
        pcg32 rng;
        pcg_generator(long long int _seed) : point_generator(_seed){
            rng(_seed);
            method_name = "pcg";
        }
        double random(){ return unif(rng); }
        void generate(vector<double> &x, region &r){
            for(unsigned i=0; i<r.dim; ++i){ 
                x[i] = r.limits[i][0] + (r.limits[i][1] - r.limits[i][0])*unif(rng); 
            }
        }
};

class mersenne_generator : public point_generator {
    public:
        std::mt19937 rng;
        mersenne_generator(long long int _seed) : point_generator(_seed){
            rng.seed(_seed);
            method_name = "mersenne";
        }
        double random(){ return unif(rng); }
        void generate(vector<double> &x, region &r){
            for(unsigned i=0; i<r.dim; ++i){ 
                x[i] = r.limits[i][0] + (r.limits[i][1] - r.limits[i][0])*unif(rng); 
            }            
        }
};

class drand48_generator : public point_generator {
    public:
        drand48_generator(long long int _seed) : point_generator(_seed){
            srand48(_seed);
            method_name = "drand48";
        }
        double random(){ return drand48(); }
        void generate(vector<double> &x, region &r){
            for(unsigned i=0; i<r.dim; ++i){ 
                x[i] = r.limits[i][0] + (r.limits[i][1] - r.limits[i][0])*drand48(); 
            }            
        }        
};

class sobol_generator : public point_generator {
    public:
        long long int sobol_seed[1] = {0};
        sobol_generator(long long int _seed) : point_generator(_seed){
            srand48(_seed);
            method_name = "sobol";
        }
        double random(){ 
            double *x;
            i8_sobol ( 1, sobol_seed, x );
            return x[0]; 
        }
        void generate(vector<double> &x, region &r){
            i8_sobol ( (int)x.size(), sobol_seed, &x[0] );
            for(unsigned i=0; i<r.dim; ++i){ 
                x[i] = r.limits[i][0] + (r.limits[i][1] - r.limits[i][0])*x[i]; 
            }            
        }            
};



#endif