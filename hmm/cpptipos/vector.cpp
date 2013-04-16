#include "tipos.h"

using namespace std;

Vector::Vector() : Matriz()  {
}

Vector::Vector(const uint m) : Matriz(1, m) {	
}

Vector::Vector(const uint m, double *d) : Matriz(1, m, d) {	
}

/*
Vector::~Vector() {
	if(a && !ref) delete[] a;
	_nRows = _nCols = 0;
}
*/

double Vector::operator() (const uint n) const {
	assert(n >= 0 && n < _nCols);

	return a[n];
}

double& Vector::operator() (const uint n) {
	//assert(m >= 0 && m < _nRows);
	assert(n >= 0 && n < _nCols);
	
	return a[n];
}

/*
Vector::Vector(const Vector &copiar) : Matriz(copiar) {
//	this->Copiar(copiar);
}

Vector& Vector::operator= (const Vector &param) {
	this->Copiar(param);
	return *(this);
}

void Vector::Copiar(const Vector &copia) {
	_nRows = copia.nRows();
	_nCols = copia.nCols();
	
	this->a = new double*[_nRows];

	for(int i=0; i<_nRows; i++) 	{
		this->a[i] = new double[_nCols];

		for(int j=0; j<_nCols; j++) {
			this->a[i][j] = copia.a[i][j];
		}
	}
}
*/

bool Vector::Guardar(const char *ruta, Vector &a) {
    ifstream archivo(ruta, ios::in);
	if(!archivo) {
		std::cerr << "Error: Archivo inexistente.";
	    return false;
    }

	int n;
	archivo >> n;
	a = Vector(n);
    
    if(n > 0) { // !archivo.failbit &&        
        for(uint i=0; i<n; ++i) {
            if(archivo.eof()) {
				std::cerr << "Error: El archivo no tiene el formato esperado.";
                return false;
            }
            archivo >> a(i);
        }
    } else {
        std::cerr << "Error: El archivo no tiene el formato esperado.";
        return false;
    }

    return true;
}
