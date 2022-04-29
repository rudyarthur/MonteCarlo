#include <iostream>
#include <vector>
#include <math.h>
#include <functional>

#include "point_generator.hpp"
#include "bins.hpp"
#include "region.hpp"
//#include "integrator.hpp"
#include "simple_mc.hpp"
#include "miser.hpp"
#include "subtraction.hpp"

using std::cout; 
using std::endl;
using std::vector;
using std::function;

double integrand(vector<double> &x){
    double sum = 1;
    //for(auto& v: x){ sum *= v;  }
    //for(auto& v: x){ sum *= exp(v);  }
    //for(auto& v: x){ sum *= sin(2*M_PI*v);  }

    /*double r = 0;
    for(auto& v: x){ r += v*v; }
    if( r < 1 ){ return 1; }
    return 0;*/

    double r = 0;
    for(auto& v: x){ r += v*v; }
    if( r < 1 && r > 0.5*0.5 ){ return 1; }
    return 0;

    return sum;
}
double exact(region &r){

    double I = 1;
    for(unsigned i=0; i<r.dim; ++i){
        //I *= ( pow(r.limits[i][1],2) - pow(r.limits[i][0],2) )*0.5;
        //I *= exp(r.limits[i][1]) - exp(r.limits[i][0]);
        //I *= -( cos(2*M_PI*r.limits[i][1]) - cos(2*M_PI*r.limits[i][0]) )/(2*M_PI);
        //I = M_PI;
        I = M_PI - M_PI/4;
    }
    return I;
}



int main(int argc, char **argv){
    

    long long int seed = 1234567;
    pcg_generator pcg(seed);

    region r( {{-1,1}, {-1,1}} );
    
    
    simple_mc mc(pcg);
    mc.integrate(integrand, r, 20000);   
    cout << "simple_mc " << mc.result_string() << " ~= " << exact(r) << endl;

    miser ms(pcg);
    ms.integrate(integrand, r, 20000);
    cout << "miser " << ms.result_string() << " ~= " << exact(r) << endl;

    dense_subtraction ds(pcg, r, 10);
    ds.integrate(integrand, r, 10000);   
    cout << "subtraction off " << ds.result_string() << " ~= " << exact(r) << endl;
    ds.initialise_approx();
    cout << "subtraction_value " << ds.subtraction_value << endl;
    ds.integrate(integrand, r, 10000);
    cout << "subtraction on " << ds.result_string() << " ~= " << exact(r) << endl;

    axis_subtraction as(pcg, r, 10);
    as.integrate(integrand, r, 10000);   
    cout << "subtraction off " << as.result_string() << " ~= " << exact(r) << endl;
    as.initialise_approx();
    cout << "subtraction_value " << as.subtraction_value << endl;
    as.integrate(integrand, r, 10000);
    cout << "subtraction on " << as.result_string() << " ~= " << exact(r) << endl;



    return 0;
}