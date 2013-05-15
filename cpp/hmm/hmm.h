#ifndef __HMM_H__
#define __HMM_H__

// #include <vector>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

typedef unsigned int uint;

class HMM {
  private: 
    static const uint DEF_MAX_ITER_ESTIM;
    static const uint DEF_MAX_ITER_HMM;

    uint MAX_ITER_ESTIM;  // max iterations allowed for EM to converge
    uint MAX_ITER_HMM;    // max iterations allowed to get best overall representant

    uint K;               // number of 'words' in dictionary
    uint T;               // number of samples in time
    uint EX;              // number of example series

    uint N;

    vector<uint> data;
    vector<uint> keys;

  public:
    class params {   
      protected: 
        matrix<double> priori;
        matrix<double> mtrans;
        matrix<double> memisn;
        
      public:
        vector<uint> observ;
        vector<uint> hidden;

        params(uint N, uint K, bool rand=false);

        friend class HMM;
        friend std::ostream& operator<< (std::ostream &out, params &p);

        static matrix<double> estimPriori(vector<uint> data); 
        static matrix<double> estimMEmisn(vector<uint> data, vector<uint> keys);
    };


    HMM(uint K, uint T, uint EX = 1, 
        uint MAX_ITER_ESTIM = HMM::DEF_MAX_ITER_ESTIM, 
        uint MAX_ITER_HMM = HMM::DEF_MAX_ITER_HMM);

    HMM(vector<uint> data,
        uint MAX_ITER_ESTIM = HMM::DEF_MAX_ITER_ESTIM, 
        uint MAX_ITER_HMM = HMM::DEF_MAX_ITER_HMM);

    void setKeys(vector<uint> keys);
    void genKeys(uint N);

    vector<uint> sample(params p, uint T);

    // static void setRNG(boost::mt19937 &rng);
    static params generateParams(uint N, uint K);
};

#endif