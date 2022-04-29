#ifndef REGION_H
#define REGION_H

#include <iostream>
#include <vector>
using std::cout; 
using std::endl;
using std::vector;

class region {
    public:

        unsigned dim;
        double volume;
        bool consistent;
        vector< vector<double> > limits;
        
        region(unsigned _dim) : dim(_dim) {
            limits.resize(dim, {0,1});
            volume = 1;
            consistent = true;
        }

        region(vector< vector<double> > _limits) : limits(_limits){ 
            dim = limits.size();
            check_consistency();
            compute_volume(); 
        }

        void check_consistency(){
            consistent = true;
            for(auto &e : limits){
                if(e[1] <= e[0]){ consistent = false; break; }
            }
        }

        inline double midpt(unsigned i){ return 0.5*(limits[i][1] + limits[i][0]); }
        inline double width(unsigned i){ return (limits[i][1] - limits[i][0]); }

        void compute_volume(){
            volume = 1;
            for(unsigned i=0; i<dim; ++i){
                volume *= width(i);
            }
        }

        void print(){
            for(unsigned i=0; i<dim; ++i){
                cout << "Dimension " << i << ": " << limits[i][0] << " " << limits[i][1] << endl;
            }
            //cout << "consistency = " << consistent << endl;
            //cout << "volume = " << volume << endl;
        }
};

#endif