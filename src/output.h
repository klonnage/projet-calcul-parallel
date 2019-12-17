#include "Vector.h"

#include <ostream>

void write_vector_to_file( const Vector &U, int Ncol, int iBegin, int iEnd, float dx, float dy, int precision = 6 );

void write_vector( const Vector &U,
                   std::ostream &outstream,
                   int           Ncol,
                   int           iBegin,
                   int           iEnd,
                   float         dx,
                   float         dy,
                   int           precision = 6 );
