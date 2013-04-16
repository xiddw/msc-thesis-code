#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <iostream>
#include <fstream>
#include <omp.h>

#include "tipos.h"

class Vector : public Matriz  {
protected:
    void Copiar(const Vector &);
public:
	Vector();
//	Vector(const Vector&);
	Vector(const uint m);
	Vector(const uint m, double *d);
	//~Vector();

	double& operator() (const uint);
	double  operator() (const uint) const;
//	Vector& operator=  (const Vector &);

	static bool Guardar(const char *ruta, Vector &);
};

#endif

