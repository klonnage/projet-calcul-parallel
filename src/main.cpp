#include "GradConj.h"
#include "Input.h"
#include "Sparse.h"
#include "Vector.h"
#include "functions.h"
#include "output.h"
#include "second_membre.h"

#include <cstring>
#include <iostream>
#include <mpi.h>
#include <algorithm>

#define str( x ) #x

using namespace std;

static double max_array(double *a, double *b, int size) {
  double max_diff = fabs(a[0] - b[0]);
  for(int i = 1; i < size; ++i) {
    max_diff = max(max_diff, fabs(a[i] - b[i]));
  }
  return max_diff;
} 

int main( int argc, char **argv )
{
    int rank, np;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &np );

    auto inputData = ReadInput( "../../data/parametres.dat" );

    int iBegin, iEnd;
    charge( inputData.rowCount,
            np,
            rank,
            iBegin,
            iEnd ); // retourne l'indice de la ligne de d�but et de celle de fin pour chaque proc


    if ( rank != 0 ) { // toujours avoir une ligne en commun pour tester schwarz
        iBegin = iBegin - inputData.coverage;
    }

    int nbLignes = iEnd - iBegin + 1;
    int nbElts   = nbLignes * inputData.colCount;

    float dX = inputData.Lrow / (inputData.rowCount+1), dY = inputData.Lcol / (inputData.colCount+1);
    int   imax = inputData.tMax / inputData.dt;

    float  a = ( ( 2. / ( dX * dX ) ) + ( 2. / ( dY * dY ) ) ) * inputData.D * inputData.dt + 1;
    float  b = (-1. / ( dX * dX )) * inputData.D * inputData.dt;
    float  c = (-1. / ( dY * dY )) * inputData.D * inputData.dt;
    Sparse A( inputData.colCount, nbLignes - 1, a, b, c );

    Vector U( nbElts ), Uprev( nbElts );
    Uprev.set_value( 0. );
    U.set_value(0.);

    // g : conditions aux bords en haut et en bas (les vecteurs conceres par MPI !!!!!),
    //   h : """""""""""""""""""" a gauche et a droite
    Vector gme( 2 * inputData.colCount );
    Vector hme( nbLignes * 2 );
    Vector termeSource( nbLignes * inputData.colCount );

    g( rank, inputData.colCount, dX, inputData.Lcol, gme, inputData.mode );
    h( rank, nbLignes, iBegin, dY, inputData.Lrow, hme, inputData.mode );
    //fsource( rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, 0, termeSource, inputData.mode );

    Vector F( nbElts );
    F.set_value( 0. );
/*    calcul_second_membre( F,
                          nbLignes,
                          inputData.colCount,
                          dX,
                          dY,
                          inputData.D,
                          inputData.dt,
                          gme,
                          hme,
                          termeSource ); // secondMembre dans second_memebre_sparse.f90*/

    MPI_Datatype line;
    MPI_Type_contiguous( inputData.colCount, MPI_DOUBLE, &line );
    MPI_Type_commit( &line );

    double *torecv = new double[inputData.colCount*inputData.coverage];
    double *comp   = new double[inputData.colCount*inputData.coverage];

    Vector _b(0);
    // boucle principale
    
    double btime = MPI_Wtime();
    for ( int i = 0; i < imax; ++i ) {
      double error;
      double error_all;

      fsource( rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, i*inputData.dt, termeSource, inputData.mode );
      calcul_second_membre( F,
                            nbLignes,
                            inputData.colCount,
                            dX,
                            dY,
                            inputData.D,
                            inputData.dt,
                            gme,
                            hme,
                            termeSource );

      do {
        _b = F;
        _b.scale( inputData.dt );
        _b += Uprev;
        GC_sparse_parallel( A, _b, U, Uprev, inputData.kmax, inputData.beta );
        Uprev = U;
        int rnext, rprev;
        if(np == 1 || np == 1) {
          break;
        }
        rnext = ( rank + 1 ) % np;
        rprev = ( rank + np - 1 ) % np;

        // First send previous last row
        double *tosend = U.data() + ( nbLignes - inputData.coverage - 1 ) * inputData.colCount;
        // Everyone send to the next rank and receive from the previous one
        MPI_Sendrecv( tosend, 1, line, rnext, 0, torecv, 1, line, rprev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        // Except 0, all update gme
        if ( rank != 0 ) { memcpy( gme.data(), torecv, inputData.colCount * sizeof( double ) ); }
        // Second row
        tosend = U.data() + (inputData.coverage)*inputData.colCount;
        MPI_Sendrecv( tosend, 1, line, rprev, 1, torecv, 1, line, rnext, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        if ( rank != np - 1 ) {
          memcpy( gme.data() + inputData.colCount, torecv, inputData.colCount * sizeof( double ) );
        }

        calcul_second_membre( F,
                              nbLignes,
                              inputData.colCount,
                              dX,
                              dY,
                              inputData.D,
                              inputData.dt,
                              gme,
                              hme,
                              termeSource );

        /* Error buffer */
        tosend = U.data() + ( nbLignes - inputData.coverage ) * inputData.colCount;
        MPI_Sendrecv( tosend, inputData.coverage, line, rnext, 2, comp, inputData.coverage, line, rprev, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        if(rank != 0) {
          error = max_array(U.data(), comp, inputData.colCount);
        } else {
          error = 0.;
        }

        MPI_Allreduce(&error, &error_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      } while ( error_all >= inputData.eps );
    }
    double diff_time = MPI_Wtime() - btime;
    double max_time;

    MPI_Reduce(&diff_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(rank == 0) {
      cout << max_time << endl;
    }

    delete[] torecv;
    delete[] comp;
    MPI_Type_free( &line );

    write_vector_to_file( U, inputData.colCount, iBegin, iEnd, dX, dY );
    MPI_Finalize();
    return 0;
}
