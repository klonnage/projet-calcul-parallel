#include "Vector.h"
#include <cblas.h>
#include <cassert>
#include <cstring>

Vector::Vector(int N) : N(N) {
    values = new double[N];
}

Vector::Vector(Vector const& v) : N(v.N), values(v.values) {

}

Vector::~Vector() {
    delete[] values;
}

double Vector::nrm2() const {
    return cblas_dnrm2(this->N, this->values, 1);
}

double Vector::dot(Vector const& v) const {
    assert(this -> N == v.N);
    return cblas_ddot(this -> N, this -> values, 1, v.values, 1);
}

double& Vector::operator[](int id) {
    return this -> values[id];
}

int Vector::size() const {
    return this -> N;
}

Vector& Vector::operator=(Vector const& v) {
    if (this -> values != nullptr) {
        delete[] this->values;
    }
    this -> N = v.N;
    this -> values = new double[this -> N];
    memcpy(this -> values, v.values, this -> N * sizeof(double));
    return *this;
}

Vector& Vector::scale(double alpha) {
    cblas_dscal(this -> N, alpha, this -> values, 1);
    return *this;
}

void axpy(double alpha, Vector const& x, Vector& y) {
    /* Authorized as x won't be modified */
    cblas_daxpy(x.size(), alpha, const_cast<Vector&>(x), 1, y, 1);
}