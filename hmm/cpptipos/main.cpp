#include <iostream>
#include <cmath>
#include "tipos.h"

using namespace std;

double fwd_back(Vector &priori, Matriz &mtrans, Matriz &memisn, Vector &data,
			  Matriz &alpha, Matriz &beta, Matriz &gamma, Matriz &xi,
			  int N, int K, int T) {

	Vector cn, ter;
	double sum, loglike;

	cn  = Vector(T);
	ter = Vector(N);

	///////////// Forward step /////////////

	for(int n=0; n<N; ++n) {
		alpha(n, 0) = priori(n);		
	}

	cn(0) = 1.0;


	for(int t=1; t<T; ++t) {		
		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {
				sum += alpha(n, t-1) * mtrans(n, m);
			}
			sum *= memisn(m, (int)data(t)-1);
			ter(m) = sum;
		}


		sum = 0.0;
		for(int n=0; n<N; ++n) {
			sum += ter(n);
		}
		cn(t) = sum;

		for(int n=0; n<N; ++n) {
			alpha(n, t) = ter(n) / cn(t);
		}
	}

	///////////// Backward step /////////////

	for(int n=0; n<N; ++n) {
		beta(n, T-1) = 1.0;
	}

	double bb;
	for(int t=T-2; t>=0; --t) {
		for(int m=0; m<N; ++m) {
			sum = 0.0;
			for(int n=0; n<N; ++n) {	
				bb = beta(n, t+1) * memisn(n, (int)data(t + 1)-1);

				xi(m, n) += mtrans(m, n) * alpha(n, t) * bb;
				sum += bb * mtrans(m, n);
			}
			beta(m, t) = sum / cn(t+1);

		}

		for(int n=0; n<N; ++n) {
			sum += beta(n, t);
		}
	}

	///////////// Log Propability /////////////
	sum = 0.0;
	for(int t=0; t<T; ++t) {
		sum += log(cn(t));
	}
	loglike = sum;

	for(int t=0; t<T; ++t) {		
		for(int n=0; n<N; ++n) {
			gamma(n, t) = alpha(n, t) * beta(n, t);
		}
	}

	return loglike;	
}

int main() {
	Matriz a, b;

	double r[] = {9, 8, 7, 6, 5, 4, 3, 2, 1};

	a = Matriz(3, 3);
	a.Llenar();
	b = Matriz(3, 3, r);

	cout << "a: " << endl << a << endl;
	cout << "b: " << endl << b << endl;
	cout << "a+b: " << endl << a+b << endl;

	return 0;
}