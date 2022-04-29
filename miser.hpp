#ifndef MISER_H
#define MISER_H


#include <limits>
#include <algorithm>

#include "integrator.hpp"
#include "region.hpp"

using std::vector;
using std::function;
using std::min;
using std::max;
double BIG = std::numeric_limits<double>::max();
double SMALL = std::numeric_limits<double>::lowest();

enum BisectionStrategy { NR, GSL };

class miser : public integrator{
    public:
        unsigned min_base_samples;  ///base case - with this many or fewer samples just do simple MC
        unsigned min_est_samples;   ///use at least this many calls to estimate variance
        double prop_samples;        ///Proportion of function calls to spend on stratification
        double dither;              ///Random noise to avoid exact bisection
        double alpha;               ///Tuning parameter, split on the basis of var^{1/(1+alpha)}
        double beta;                ///beta = 1/(1+alpha)
        int last_split;             ///remember the last co=ordinate axis that was split
        double var_threshold;       ///threshold for no variance (for constant functions)
        BisectionStrategy strategy; ///how to decide on bisection

        miser(point_generator &_p) : integrator(_p) {
            min_base_samples = 60;
            min_est_samples = 15;
            prop_samples = 0.1;
            dither = 0;
            alpha = 2;
            beta = 2./(1+alpha);
            last_split = -1;
            var_threshold = 1e-10;
            strategy = NR;
        }

        /*This is the strategy used in Numerical Recipies 
            Variance of f in a region = max(f) - min(f)
        */
        inline vector< vector<double> > bisection_minmax(function<double(vector<double>&)> f, unsigned Np, region &r, vector<double> &mid){

            vector< vector<double> > bi_variances(r.dim, {BIG,SMALL,BIG,SMALL});
            vector<double> x(r.dim);
            for(unsigned i=0; i<Np; ++i) { //compute the variances on either side of bisection
                p.generate(x, r);
                double val = f(x); 
                ++nevals;
                for(unsigned i=0; i<x.size(); ++i){
                    if(x[i] < mid[i]){
                        bi_variances[i][0] = min(bi_variances[i][0], val);
                        bi_variances[i][1] = max(bi_variances[i][1], val);
                    } else {
                        bi_variances[i][2] = min(bi_variances[i][2], val);
                        bi_variances[i][3] = max(bi_variances[i][3], val);
                    }
                }
            }

            return bi_variances;

        }

        /*This is the strategy used by GSL
            Variance of f in a region = <f^2> - <f>^2
        */
        inline vector< vector<double> > bisection_variance(function<double(vector<double>&)> f, unsigned Np, region &r, vector<double> &mid){

            dim = r.dim;
            
            vector<bin> left(r.dim);
            vector<bin> right(r.dim);
            vector<double> x(r.dim);
            for(unsigned i=0; i<Np; ++i) { //compute the variances on either side of bisection
                p.generate(x, r);
                double val = f(x); 
                ++nevals;
                for(unsigned i=0; i<x.size(); ++i){
                    if(x[i] < mid[i]){
                        left[i].accumulate(val);
                    } else {
                        right[i].accumulate(val);
                    }
                }
            }
            
            vector< vector<double> > bi_variances(r.dim, {0,0});
            for(unsigned i=0; i<x.size(); ++i){
                left[i].normalise();
                right[i].normalise();

                bi_variances[i][0] = sqrt( (left[i].sum2 - left[i].sum*left[i].sum)/left[i].count )*(mid[i] - r.limits[i][0])/r.width(i);
                bi_variances[i][1] = sqrt( (right[i].sum2 - right[i].sum*right[i].sum)/right[i].count )*(r.limits[i][1] - mid[i])/r.width(i);

            }
            return bi_variances;

        }



        integrate_result integrate(function<double(vector<double>&)> f, region &r, unsigned N){

            vector<double> x(r.dim); 

            //base case
            if(N < min_base_samples){

                return accumulate(f, r, N);
                
            } else {
                vector<double> mid(r.dim);
                //bisect the axes
                for(unsigned i=0; i<r.dim; ++i){ mid[i] = r.midpt(i) + dither*(-1 + 2*p.random())*r.width(i)*0.5;  }
                

                unsigned Np = max( (unsigned)(N*prop_samples), min_est_samples );
                vector< vector<double> > bi_variances; 
                switch(strategy){
                    case NR: bi_variances = bisection_minmax(f, Np, r, mid); break;
                    case GSL: bi_variances = bisection_variance(f, Np, r, mid); break;
                }
                N -= Np;

                //now find the dimension where the variance sum is smallest
                //NR and GSL don't seem to account for ties! So gather all the results and check at the end.
                int split_dim = -1;
                vector<double> dim_variance(r.dim); 
                vector< vector<double> > lr_variance(r.dim, {0,0});
                for(unsigned i=0; i<r.dim; ++i) {
                    double s_left, s_right;
                    switch(strategy){
                        case NR:
                            s_left  = pow(bi_variances[i][1] - bi_variances[i][0], beta);
                            s_right = pow(bi_variances[i][3] - bi_variances[i][2], beta);
                            break;
                        case GSL:
                            s_left  = pow(bi_variances[i][0], beta);
                            s_right = pow(bi_variances[i][1], beta);
                            break;
                    }
                    dim_variance[i] = s_left + s_right;
                    lr_variance[i][0] = s_left;
                    lr_variance[i][1] = s_right;
                }
                double min_variance = *std::min_element(dim_variance.begin(), dim_variance.end());
                //if(min_variance < var_threshold){ return accumulate(f, r, N); } //very small variance region, bail here?

                vector<double> min_var; for(unsigned i=0; i<r.dim; ++i){ if(dim_variance[i] == min_variance){ min_var.push_back(i); } }
                if(min_var.size() == 1){  
                    split_dim = min_var[0]; 
                } else{
                    int idx = (int)(min_var.size() * p.random());           
                    split_dim = min_var[ idx ];
                    //avoid splitting the same axis twice in a row if you can help it...
                    while(last_split == split_dim){split_dim = min_var[ (int)(min_var.size() * p.random()) ]; } 
                }
                last_split = split_dim; 
                

                //calculate the number of points per region
                //subvolume size AND variance determine point split                
                //size of the left subvolume always 1/2 if dither = 0
                double fraction_left = (mid[split_dim] - r.limits[split_dim][0])/r.width(split_dim); 
                double fraction_right = 1 - fraction_left;
                double sigma_left = lr_variance[split_dim][0];
                double sigma_right = lr_variance[split_dim][1];
                double left_ratio;
                double right_ratio;
                if(sigma_left == 0 && sigma_right == 0){ //this shouldn't happen!
                    left_ratio = fraction_left / (fraction_right + fraction_left);
                    right_ratio = fraction_right / (fraction_right + fraction_left);
                } else {
                    sigma_left *= fraction_left;
                    sigma_right *= fraction_right;
                
                    left_ratio = sigma_left / (sigma_left + sigma_right);
                    right_ratio = sigma_right / (sigma_left + sigma_right);
                }
                unsigned N_left = min_est_samples + (N - 2*min_est_samples)*left_ratio;  
                unsigned N_right = min_est_samples + (N - 2*min_est_samples)*right_ratio; 

                //make new regions and recurse
                region r_left(r); r_left.limits[split_dim][1] = mid[split_dim]; r_left.compute_volume();
                integrate_result res_left = integrate(f, r_left, N_left);

                region r_right(r); r_right.limits[split_dim][0] = mid[split_dim]; r_right.compute_volume();
                integrate_result res_right = integrate(f, r_right, N_right);

                //The Var(<f>') = (1/4)*[ Var_a(f)/Na + Var_b(f)/Nb] formula from Numerical Recipies
                res.average = fraction_left*res_left.average + fraction_right*res_right.average;
                res.integral = res_left.integral + res_right.integral;
                res.variance = fraction_left*fraction_left*res_left.variance + fraction_right*fraction_right*res_right.variance;
                res.error =  sqrt( res_left.error*res_left.error + res_right.error*res_right.error) ;

                return res;
            }  
        }
};

#endif