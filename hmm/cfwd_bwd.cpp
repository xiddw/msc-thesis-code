#include <stdlib.h>
#include <math.h>

#include "mex.h"
#include "cpptipos/tipos.h"

double fwd_back(Vector &priori, Matriz &mtrans, Matriz &memisn, Vector &data,
			  Matriz &alpha, Matriz &beta, Matriz &gamma, Matriz &xi, Vector &cn,
			  int N, int K, int T) {

	double sum, loglike;
	double aa, bb;

	// Vector cn  = Vector(T);
	Matriz xx  = Matriz(N, N);

	///////////// Forward step /////////////
	for(int n=0; n<N; ++n) {
		alpha(n, 0) = priori(n) * memisn(n, (int)data(0)-1);
	}
	cn(0) = alpha.colNorm(0);

	for(int t=1; t<T; ++t) {		
		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {
				aa = alpha(n, t-1) * memisn(m, (int)data(t)-1);			
				sum += aa * mtrans(n, m);
			}
			alpha(m, t) = sum;
		}

		cn(t) = alpha.colNorm(t);
	}

	///////////// Backward step /////////////
	for(int n=0; n<N; ++n) {
		beta(n, T-1) = 1.0;
	}
	
	for(int t=T-2; t>=0; --t) {
		xx = 0.0;

		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {	
				bb = beta(n, t+1) * memisn(n, (int)data(t + 1)-1);			
				sum += bb * mtrans(m, n);

				xx(m, n) += mtrans(m, n) * alpha(m, t) * bb;
			}
			beta(m, t) = sum;
		}
		xx.Norm();
		xi = xi + xx;

		beta.colNorm(t);
	}

	xi.Norm();

	///////////// Log Propability /////////////
	sum = 0.0;
	// mexPrintf("AAAAAA: \n");
	for(int t=0; t<T; ++t) {
		sum += log(cn(t));
		// mexPrintf("%f, ", cn(t));
	}
	loglike = sum;
	// mexPrintf("\n:BBBBBB \n");

	///////////// Gamma calc /////////////
	for(int t=0; t<T; ++t) {		
		for(int n=0; n<N; ++n) {
			gamma(n, t) = alpha(n, t) * beta(n, t);
		}

		gamma.colNorm(t);
	}

	return loglike;	
}

void mexFunction(int nlhs, mxArray* plhs[], 			// Variables de entrada
				 int nrhs, const mxArray *prhs[]) {		// Variables de salida

	if(nrhs != 4) {
		mexErrMsgTxt("Se deben recibir 4 parametros de entrada");
	}

	if(nlhs != 5+1) {
		mexErrMsgTxt("Se deben recibir 5 parametros de salida");
	}

	const int 
		N  = (int)mxGetM(prhs[2]),	// Numero de speakers
		K  = (int)mxGetN(prhs[2]),	// Numero de palabras en el diccionario
		EX = (int)mxGetM(prhs[3]),	// Numero de ejemplos
		T  = (int)mxGetN(prhs[3]);	// Numero de muestras en el tiempo

	Vector priori = Vector(N, mxGetPr(prhs[0]));
	Matriz mtrans = Matriz(N, N, mxGetPr(prhs[1]));
	Matriz memisn = Matriz(N, K, mxGetPr(prhs[2]));
	Vector data   = Vector(T, mxGetPr(prhs[3]));

	plhs[0] = mxCreateDoubleMatrix(N, T, mxREAL);	// alpha
	plhs[1] = mxCreateDoubleMatrix(N, T, mxREAL);	// beta
	plhs[2] = mxCreateDoubleMatrix(N, T, mxREAL);	// gamma
	plhs[3] = mxCreateDoubleMatrix(N, N, mxREAL);	// xi
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);	// loglike
	
	plhs[5] = mxCreateDoubleMatrix(T, 1, mxREAL);	// cn (temp)

	Matriz alpha = Matriz(N, T, mxGetPr(plhs[0]));
	Matriz beta  = Matriz(N, T, mxGetPr(plhs[1]));
	Matriz gamma = Matriz(N, T, mxGetPr(plhs[2]));
	Matriz xi    = Matriz(N, N, mxGetPr(plhs[3]));
	
	Vector cn    = Vector(T, mxGetPr(plhs[5]));
	
	double loglike = fwd_back(priori, mtrans, memisn, data,
			 alpha, beta, gamma, xi, cn,  
			 N, K, T);

	*mxGetPr(plhs[4]) = loglike;
}