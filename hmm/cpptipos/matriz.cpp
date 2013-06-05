#include <iostream>
#include <cstdlib>

#include "tipos.h"

using namespace std;

uint Matriz::nCols() const { return _nCols; }
uint Matriz::nRows() const { return _nRows; }

Matriz::Matriz() {
	_nRows = _nCols = 0;
	ref = false;
}

Matriz::Matriz(const uint m) {
	ref = false;
	this->Inicializar(m, m);
}

Matriz::Matriz(const uint m, const uint n) {
	ref = false;
	this->Inicializar(m, n);
}

Matriz::Matriz(const uint m, const uint n, double *d) {
	ref = true;

	_nRows = m; 
	_nCols = n;

	a = d;
}

Matriz::Matriz(const Matriz &copiar) {
	ref = false;
	this->Copiar(copiar);
}


Matriz::~Matriz() {
	if(a && !ref) delete[] a;
	_nRows = _nCols = 0;
}	

ostream& operator<<(ostream &os, const Matriz &m) {	
	int r = m.nRows();
	int c = m.nCols();

	for(uint i=0; i<r; ++i) {
		for(uint j=0; j<c; ++j) {
			os.precision(6);
			os.width(12);
			// os << m(i, j) << " ";
			// os << m.a[m.nRows()*j + i];
			os << m.a[m.nCols()*i + j];
		}
		os << std::endl;
	}
	os << std::endl;
	
	return os;
}

Matriz operator+ (const Matriz &a, const Matriz &b) {
	uint r = a.nRows();
	uint c = a.nCols();

	Matriz temp = Matriz(r, c);

	for(uint i=0; i<r; ++i) {
		for(uint j=0; j<c; ++j) {
			temp(i, j) = a(i, j) + b(i, j);
		}
	}

	return (temp);
}

Matriz operator- (const Matriz &a, const Matriz &b) {
	uint r = a.nRows();
	uint c = a.nCols();

	Matriz temp = Matriz(r, c);

	for(uint i=0; i<r; ++i) {
		for(uint j=0; j<c; ++j) {
			temp(i, j) = a(i, j) - b(i, j);
		}
	}

	return (temp);
}

Matriz& Matriz::operator= (const Matriz &param) {
	this->Copiar(param);
	return *(this);
}

Matriz& Matriz::operator= (const double v) {
	this->Llenar(v);
	return *(this);
}

inline double Matriz::operator() (const uint m, const uint n) const {
	assert(m >= 0 && m < _nRows);
	assert(n >= 0 && n < _nCols);
	
	return a[_nRows*n + m];
}

inline double& Matriz::operator() (const uint m, const uint n) {
	assert(m >= 0 && m < _nRows);
	assert(n >= 0 && n < _nCols);
	
	return a[_nRows*n + m];
}

bool Matriz::EsCuadrada() const {
	return (_nCols == _nRows);
}

void Matriz::Llenar() {
	for(uint i=0; i<_nRows; ++i) {
		for(uint j=0; j<_nCols; ++j) {
			(*this)(i, j) = (rand()%20 + 1);
		}
	}
}

void Matriz::Llenar(const double v) {
	for(uint i=0; i<_nRows; ++i) {
		for(uint j=0; j<_nCols; ++j) {
			(*this)(i, j) = v;
		}
	}
}

bool Matriz::Guardar(const char *ruta, Matriz &a) {
    ifstream archivo(ruta, ios::in);
	if(!archivo) {
		std::cerr << "Error: Archivo inexistente.";
	    return false;
    }

	int m, n;
	archivo >> m;
	archivo >> n;
	a = Matriz(m, n);
    
    if(m > 0) { // !archivo.failbit &&        
        for(uint i=0; i<m; ++i) {
            for(uint j=0; j<n; ++j) {
                if(archivo.eof()) {
					std::cerr << "Error: El archivo no tiene el formato esperado.";
                    return false;
                }
                archivo >> a(i, j);
            }
        }
    } else {
        std::cerr << "Error: El archivo no tiene el formato esperado.";
        return false;
    }

    return true;
}


/* 
 * Métodos privados 
 */

void Matriz::Inicializar(const uint m, const uint n) {
	_nRows = m;
	_nCols = n;
	
    try {
        this->a = new double[_nRows * _nCols];
    } catch(std::bad_alloc& exc) {
        return;
    }

	this->Llenar(0.0);
}

void Matriz::Copiar(const Matriz &copia) {
	try {
    	if(_nRows != copia.nRows() || 
    	   _nCols != copia.nCols()) {

	    	_nRows = copia.nRows();
			_nCols = copia.nCols();

	        this->a = new double[_nRows * _nCols];
	    }
    } catch(std::bad_alloc& exc) {
        return;
    }
    
	for(uint i=0; i<_nRows; ++i) { 
		for(uint j=0; j<_nCols; ++j) {
			(*this)(i, j) = copia(i, j);
		}
	}
}

double Matriz::rowNorm(const uint m) {
	assert(m >= 0 && m < _nRows);

	double s = 0.0;
	for(uint j=0; j<_nCols; ++j) {
		s += (*this)(m, j);
	}

	s += (s == 0.0);

	for(uint j=0; j<_nCols; ++j) {
		(*this)(m, j) /= s;
	}

	return s;
}

double* Matriz::rowNorm(void) {
	double *scale = new double[_nRows];

	for(uint i=0; i<_nRows; ++i) {
		scale[i] = this->rowNorm(i);
	}

	return scale;
}

double Matriz::colNorm(const uint n) {	
	assert(n >= 0 && n < _nCols);

	double s = 0.0;
	for(uint i=0; i<_nRows; ++i) {
		s += (*this)(i, n);
	}

	s += (s == 0.0);

	for(uint i=0; i<_nRows; ++i) {
		(*this)(i, n) /= s;
	}

	return s;
}

double* Matriz::colNorm(void) {
	double *scale = new double[_nCols];

	for(uint j=0; j<_nCols; ++j) {
		scale[j] = this->colNorm(j);
	}

	return scale;
}

double Matriz::Norm(void) {
	double s = 0.0;

	for(uint i=0; i<_nRows; ++i) {
		for(uint j=0; j<_nCols; ++j) {
			s += (*this)(i, j);
		}
	}

	s += (s == 0.0);

	for(uint i=0; i<_nRows; ++i) {
		for(uint j=0; j<_nCols; ++j) {
			(*this)(i, j) /= s;
		}
	}

	return s;
}