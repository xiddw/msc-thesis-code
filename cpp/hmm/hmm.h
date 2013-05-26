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
  public:
    class params {   
      protected: 
        matrix<double> priori;
        matrix<double> mtrans;
        matrix<double> memisn;

        uint N;           // number of hidden states
        uint K;           // number of 'words' in dictionary        
        
      public:
        vector<uint> observ;
        vector<uint> hidden;

        params(uint N, uint K, bool rand=false);

        friend class HMM;
        friend std::ostream& operator<< (std::ostream &out, params &p);

        static matrix<double> estimPriori(vector<uint> observ); 
        static matrix<double> estimMEmisn(vector<uint> observ, vector<uint> hidden);
        static vector<uint> estimHidden(vector<uint> observ, vector<uint> keys);
    };

  private: 
    static const uint DEF_MAX_ITER_ESTIM;
    static const uint DEF_MAX_ITER_HMM;

    uint MAX_ITER_ESTIM;  // max iterations allowed for EM to converge
    uint MAX_ITER_HMM;    // max iterations allowed to get best overall representant

    uint K;               // number of 'words' in dictionary
    uint T;               // number of samples in time
    uint EX;              // number of example series
    uint N;               // number of hidden states

    vector<uint> data;
    // vector<uint> keys;

    bool ConvergedEM(double ll1, double ll0, double threshold = 1e-10, bool HasIncresed = true);
    double CalculateValues(params in, vector<uint> data, params &out);
    double BackwardForward(params in, vector<uint> data, 
                           matrix<double> &alpha, matrix <double> &beta,
                           matrix<double> &gamma, matrix <double> &xi);


  public:
    HMM(uint K, uint T, uint EX = 1, 
        uint MAX_ITER_ESTIM = HMM::DEF_MAX_ITER_ESTIM, 
        uint MAX_ITER_HMM = HMM::DEF_MAX_ITER_HMM);

    HMM(vector<uint> data,
        uint MAX_ITER_ESTIM = HMM::DEF_MAX_ITER_ESTIM, 
        uint MAX_ITER_HMM = HMM::DEF_MAX_ITER_HMM);

    vector<uint> sample(params p, uint T);

    void EM(params &p, uint MAX_ITER_HMM = 0);

    static vector<uint> generateKeys(uint N, uint K);
    static params generateParams(uint N, uint K);
};

#endif