#pragma once

class Vector {
    double *values;
    int N;
public:
    /* Par défaut : vecteur null */
    Vector() : N(0), values(nullptr) {}
    Vector(int N);

    /* Copies superficielles */
    Vector(Vector const& v);
    void set(Vector const& v);

    ~Vector();
    /* Copie profonde de v */
    Vector& operator=(Vector const& v);

    Vector& operator+=(Vector const& v);
    Vector& operator-=(Vector const& v);

    Vector& scale(double alpha);

    double& operator[](int id);
    double  dot(Vector const& v)  const;
    double  nrm2() const;

    int size() const;
    operator double*() {return values;};
    operator double const*() {return const_cast<double const*>(values);};
};

Vector operator+(Vector const &u, Vector const &v);
Vector operator-(Vector const &u, Vector const &v);
void axpy(double alpha, Vector const& x, Vector& y);