#include "Vector.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <mkl_cblas.h>

Vector::Vector( int N )
    : N( N )
{
    values = new double[N];
}

Vector::Vector( Vector const &v )
    : N( v.N )
    , values( v.values )
{
}

Vector::~Vector()
{
    delete[] values;
    values = nullptr;
}

Vector &Vector::operator+=( Vector const &v )
{
    cblas_daxpy( std::min( size(), v.size() ), 1., v.data(), 1, values, 1 );
    return *this;
}

double Vector::nrm2() const
{
    return cblas_dnrm2( this->N, this->values, 1 );
}

double Vector::dot( Vector const &v ) const
{
    assert( this->N == v.N );
    return cblas_ddot( this->N, this->values, 1, v.values, 1 );
}

int Vector::size() const
{
    return this->N;
}

Vector &Vector::operator=( Vector const &v )
{
    if ( this->values != nullptr && N != 0 ) { delete[] this->values; }
    this->N      = v.N;
    this->values = new double[this->N];
    memcpy( this->values, v.values, this->N * sizeof( double ) );
    return *this;
}

Vector &Vector::operator-=( Vector const &v )
{
    cblas_daxpy( std::min( size(), v.size() ), -1., v.data(), 1, values, 1 );
    return *this;
}

Vector &Vector::scale( double alpha )
{
    cblas_dscal( this->N, alpha, this->values, 1 );
    return *this;
}

void Vector::set_value( double v )
{
    if ( v == 0. ) { memset( this->values, 0, sizeof( double )*N ); }
    else {
        for ( int i = 0; i < this->N; ++i ) { this->values[i] = v; }
    }
}

void Vector::set( Vector const &v ) {
    N = v.N;
    memcpy(values, v.values, N);
}

void axpy( double alpha, Vector const &x, Vector &y )
{
    /* Authorized as x won't be modified */
    cblas_daxpy( x.size(), alpha, x.data(), 1, y.data(), 1 );
}

std::ostream& operator<<(std::ostream& os, Vector const& v) {
  os << "{ ";
  for(int i = 0; i < v.size(); ++i) {
    os << v[i] << " "; 
  }
  os << "}";
  return os;
}
