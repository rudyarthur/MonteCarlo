#ifndef BINS_H
#define BINS_H

#include <vector>
#include <iostream>
using std::cout; 
using std::cerr; 
using std::endl;
using std::vector;

class bin{
    public:
        unsigned count;
        double sum;
        double sum2;

        bin() {
            count = 0;
            sum = 0;
            sum2 = 0;
        }        

        void accumulate(double val){
            ++count;
            sum += val;
            sum2 += val*val;
        }

        void normalise(){
            sum /= count;
            sum2 /= count;
        }

        void print(){
            cout << count << " " << sum << " " << sum2 << endl;
        }
};

class bin_axis{
    public:
        vector<double> edge;
        double a, b;
        unsigned B;

        bin_axis(){a = 0; b = 0; B = 0; edge = {};};

        bin_axis(unsigned _B, double _a=0, double _b=1) : B(_B), a(_a), b(_b){        
            make_uniform_bins();
        }

        void make_uniform_bins(){
            edge.resize(B+1, 0);
            edge[0] = a;
            for(unsigned i=1; i<B; ++i){
                edge[i] = edge[i-1] + (b-a)/B;
            } edge[B] = b;
        }

        inline double width(unsigned i){ return edge[i+1] - edge[i]; }

        int bin_id(double x){
            if( x < a || x > b){ return -1; }
            for(int i=1; i<B+1; ++i){ 
                if(x < edge[i]){return i-1;} 
            }
            cerr << "logic bomb!" << endl;
            exit(1);
        }

        void print(){
            for(int i=0; i<B; ++i){
                cout << "bin " << i << " = [ " << edge[i] << " , " << edge[i+1] << "]" << endl;
            }
        }
};

/*
class bin_axis {
    public:
        unsigned NB;
        vector<unsigned> count;
        vector<double> sum;
        vector<double> sum2;
        vector<double> edges;
        vector<double> widths;
        double total_sum;
        unsigned total_count;

        double a, b;
        bin_axis(unsigned _NB, double _a=0, double _b=1) : NB(_NB), a(_a), b(_b){        
            clear();
            uniform_bins();
        }
        void clear(){
            count.resize(NB, 0);
            sum.resize(NB, 0);
            sum2.resize(NB, 0); 
            total_sum = 0;
            total_count = 0;
        }
        void uniform_bins(){
            edges.resize(NB+1, 0);
            widths.resize(NB, (b-a)/NB);
            edges[0] = a;
            for(unsigned i=1; i<NB; ++i){
                edges[i] = edges[i-1] + (b-a)/NB;
            } edges[NB] = b;
        }
        
        void accumulate(double x, double val){
            unsigned bid = bin_id(x);
            ++count[bid];
            sum[bid] += val;
            sum2[bid] += val*val;
            ++total_count;
            total_sum += val;
        }
        void print(){
            for(int i=0; i<NB; ++i){
                cout << count[i] << " " << sum[i] << " " << sum2[i] << endl;
            }
        }
        int bin_id(double x){
            if( x < a || x > b){ return -1; }
            for(int i=1; i<NB+1; ++i){ 
                if(x < edges[i]){return i-1;} 
            }
            cerr << "logic bomb!" << endl;
            exit(1);
        }

       void normalise(){
           for(unsigned i=0; i<NB; ++i){
               sum[i] /= count[i];
               sum2[i] = (sum2[i]/count[i]) - (sum[i]*sum[i]);
           }
           total_sum /= total_count;
       }
        
};*/

#endif