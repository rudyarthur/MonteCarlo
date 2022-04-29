#ifndef SIMPLE_MC_H
#define SIMPLE_MC_H

#include <vector>
#include <functional>

#include "integrator.hpp"

using std::vector;
using std::function;

class simple_mc : public integrator{
    public:
        simple_mc(point_generator &_p) : integrator(_p) {}
        integrate_result integrate(function<double(vector<double>&)> f, region &r, unsigned N){ return accumulate(f, r, N); }  
};

#endif