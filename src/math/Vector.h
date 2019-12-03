#pragma once

class Vector {
    double *values;
    int N;
public:
    Vector() : N(0), values(nullptr);
    Vector(int N);
    Vector(Vector const& v);
    ~Vector();
    Vector& operator=(Vector const& v);

    Vector& operator*(double alpha);
    Vector& operator+=(Vector const& v);
    Vector& operator-=(Vector const& v);

    double& operator[](int id);
    double  dot()  const;
    double  nrm2() const;

    int size() const;
    operator double*() {return values;};
    operator double const*() {return const_cast<double const*>(values);};
};

Vector& operator+(Vector const &u, Vector const &v);
Vector& operator-(Vector const &u, Vector const &v);
void axpy(double a, Vector const* x, Vector* y);