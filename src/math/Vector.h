#pragma once
#include <iostream>

class Vector {
  int     N;
  double *values;

  public:
    /* Par défaut : vecteur null */
    Vector()
        : N( 0 )
        , values( nullptr )
    {
    }
    Vector( int N );

    /* Copies superficielles */
    Vector( Vector const &v );
    void set( Vector const &v );

    void set_value( double v );

    ~Vector();
    /* Copie profonde de v */
    Vector &operator=( Vector const &v );

    Vector &operator+=( Vector const &v );
    Vector &operator-=( Vector const &v );

    Vector &scale( double alpha );

    double *      data() { return values; }
    const double *data() const { return values; }
    double &      operator[]( int id ) { return values[id]; }

    double        operator[]( int id ) const { return values[id]; };
    double        dot( Vector const &v ) const;
    double        nrm2() const;

    int size() const;
};

std::ostream& operator<<(std::ostream& os, Vector const& v);
Vector operator+( Vector const &u, Vector const &v );
Vector operator-( Vector const &u, Vector const &v );
void   axpy( double alpha, Vector const &x, Vector &y );
