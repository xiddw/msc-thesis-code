#ifndef __MATRIZ_H__
#define __MATRIZ_H__

#include <ostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <omp.h>

#include "tipos.h"

class Matriz {
protected: 
	int _nCols;
	int _nRows;

	double *a;
	void Inicializar(const uint m, const uint n);
	void Copiar(const Matriz &);

	bool ref;

	std::string _title;

public:
	uint nCols() const;
	uint nRows() const;

	Matriz();
	Matriz(const uint m);
	Matriz(const uint m, const uint n);
	Matriz(const uint m, const uint n, double *d);

	Matriz(const Matriz&);	
	~Matriz();

	friend Matriz operator+ (const Matriz &, const Matriz &);
	friend Matriz operator- (const Matriz &, const Matriz &);
	friend std::ostream& operator<<(std::ostream &os, const Matriz &);

	Matriz& operator=  (const double v);
	Matriz& operator=  (const Matriz &);
	double& operator() (const uint, const uint);
	double	operator() (const uint, const uint) const;

	bool EsCuadrada() const;
	void Llenar();
	void Llenar(const double v);

	double rowNorm(const uint m);
	double colNorm(const uint n);

	double Norm(void);

	double* rowNorm(void);
	double* colNorm(void);

	static bool Guardar(const char *ruta, Matriz &);
};

#endif