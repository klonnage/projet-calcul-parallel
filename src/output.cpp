#include "output.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <mpi.h>
#include <sstream>
#include <string>

void write_vector_to_file( const Vector &U, int Ncol, int iBegin, int iEnd, float dx, float dy, int precision /*= 6 */ )
{
    // get file name
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    std::time_t       writeTime = std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() );
    std::stringstream ss;
    ss << std::ctime( &writeTime );
    std::string fileName( "out/GradConj_" );
    fileName += std::to_string( rank );
    fileName += '_';
    fileName += ss.str();
    fileName += std::string( ".txt" );
    auto isnon_alphanumeric = []( char c ) { return std::isspace( c ) || std::isalnum( c ); };
    std::replace_if( fileName.begin(), fileName.end(), isnon_alphanumeric, '_' );

    // create file
    std::ofstream outputFile( fileName, std::ios::out | std::ios::trunc );
    if ( !outputFile.is_open() ) {
        std::cerr << "[Rank " << rank << "] could not create file " << fileName
                  << ".\n Maybe no \"out\" directory exists in the same directory as exe ?" << std::endl;
        return;
    }

    // write to file
    write_vector( U, outputFile, Ncol, iBegin, iEnd, dx, dy, precision );

    // clean
    outputFile.close();

    return;
}

void write_vector( const Vector &U,
                   std::ostream &outstream,
                   int           Ncol,
                   int           iBegin,
                   int           iEnd,
                   float         dx,
                   float         dy,
                   int           precision /*= 6*/ )
{
    auto oldPrecision = outstream.precision( precision );
    for ( int row = iBegin; row < iEnd; ++row ) {
        float x = row * dx;
        for ( int j = 0; j < Ncol; ++j ) {
            float y = j * dy;
            outstream << x << ' ' << y << ' ' << U[row * Ncol + j] << '\n';
        }
    }
    outstream.precision( precision );
}