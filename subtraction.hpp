#ifndef SUBTRACTION_H
#define SUBTRACTION_H

#include <vector>
#include <functional>

#include "integrator.hpp"

using std::vector;
using std::function;

class subtraction : public integrator {
    public:
        vector<bin_axis> ba;
        unsigned B;
        double subtraction_value;
        bool subtract;

        subtraction(point_generator &_p, region &r, unsigned _B) : integrator(_p), B(_B) {
            dim = r.dim;
            ba.resize(0);
            for(unsigned i=0; i<dim; i++){
                ba.push_back( bin_axis(B, r.limits[i][0], r.limits[i][1]) );
            }
            subtraction_value = 0;
            subtract = false;
        }


        integrate_result integrate(function<double(vector<double>&)> f, region &r, unsigned N){
                double sum = 0;
                double sum2 = 0;
                vector<double> x(r.dim); 

                for(unsigned i=0; i<N; ++i){
                    p.generate(x, r);
                    double val = f(x); 
                    val -= (subtract) ? approx(x) : 0;
                    ++nevals;
                    sum += val;
                    sum2 += val*val;
                    accumulate_in_bins(x, val);
                }

                res.normalise(sum, sum2, r.volume, N);
                res.integral += subtraction_value;

                return res;
            }  

            virtual void initialise_approx() = 0;
            virtual double approx(vector<double> &x) = 0;
            virtual void accumulate_in_bins(vector<double> &x, double val) = 0;
};

class axis_subtraction : public subtraction {
    public:

        vector< vector<bin> > bins;
        vector< vector<bin> > sub;

        axis_subtraction(point_generator &_p, region &r, unsigned _B) : subtraction(_p, r, _B) {
            bins.resize( r.dim, vector<bin>(B) );      
            sub.resize( r.dim, vector<bin>(B) );      
        }

        void accumulate_in_bins(vector<double> &x, double val){
            for(unsigned j=0; j<x.size(); ++j){ bins[j][ ba[j].bin_id(x[j]) ].accumulate(val); }
        }

        ///Exact value of approx integral
        //I = int_0^1 f(x) g(y) h(z) = Ix Iy Iz = <f> <g> <h> |dx| |dy| |dz| 
        ///sx_i ~ <f>_i <g> <h>
        ///sum_i dx_i sx_i ~ <g> <h> Ix
        ///sum_i dx_i sx_i/sum_i dx_i ~ <g> <h> <f>
        ///(sum_i  dx_i sx_i) * (sum_j dy_j sy_j) * (sum_j dz_k sz_k) = Iz <g> <h> Iy <f> <h> Iz <f> <g> 
        /// = I (I/V)^{dim-1}
        void initialise_approx(){ 
            sub = bins;
            for(auto &axis : sub){ for(auto &e : axis){ e.normalise(); } }

            subtraction_value = 1;
            for(unsigned d=0; d<dim; ++d){
                double axis_sum = 0;
                for(unsigned i=0; i<B; ++i){axis_sum += sub[d][i].sum * ba[d].width(i); }
                subtraction_value *= axis_sum;
            }
            subtraction_value /= pow(res.average,dim-1);
            subtract=true;
        }

        //TODO volume factors!
        //I = int_0^1 f(x) g(y) h(z) = Ix Iy Iz
        //sx_i * sy_j * sz_k = 
        //( <f>_i Iy Iz ) * ( <g>_j Ix Iz ) * ( <h>_k Ix Iy ) =
        //=> f(r) ~ <f>_i <g>_j <h>_k ~ ( prod_i sd_i ) / I^{dim-1}
        double approx(vector<double> &x){
            double out = 1;
            for(unsigned i=0; i<dim; ++i){ 
                out *= sub[i][ ba[i].bin_id(x[i]) ].sum; 
            }
            return out/pow(res.average,dim-1);
        }

        
};

class dense_subtraction : public subtraction {
    public:
        unsigned NB;
        vector<bin> bins;

        dense_subtraction(point_generator &_p, region &r, unsigned _B) : subtraction(_p, r, _B) {
            NB = pow(B, r.dim);
            bins.resize(NB);
        }

        ///indexing function -> find the bin index containing the point x
        ///[7,2] = 7 + 10*2 = 27
        inline unsigned get_bin_idx(vector<double> &x){
            unsigned idx = 0;

            unsigned tens = 1;
            for(unsigned d=0; d<x.size(); ++d){ 
                idx += ba[d].bin_id(x[d]) * tens;
                tens *= B; 
            }

            return idx;
        }

        ///indexing function -> find the bin coords correcponding to the index idx
        ///27 = [7,2]
        inline vector<unsigned> idx_to_bin(unsigned idx){
            vector<unsigned> bdim(ba.size());

            for(unsigned d=0; d<bdim.size(); ++d){ 
                bdim[d] = idx % B; 
                idx /= B; 
            }

            return bdim;
        }

        void accumulate_in_bins(vector<double> &x, double val){ 
            bins[ get_bin_idx(x) ].accumulate(val); 
        }

        ///Exact value of approximate integral
        ///sum_{b in bins} vol(b) <f>_b
        void initialise_approx(){ 
            subtraction_value = 0;

            for(unsigned i=0; i<NB; ++i){
                bins[i].normalise(); 
                vector<unsigned> bdim = idx_to_bin(i);
                double vol = 1;
                for(unsigned d=0; d<bdim.size(); ++d){ vol *= ba[d].width( bdim[d] ); }
                subtraction_value += vol*bins[i].sum;
            } 

            subtract=true;
        }

        ///Function approximation in bin containing x
        double approx(vector<double> &x){
            return bins[ get_bin_idx(x) ].sum;
        }
};



#endif