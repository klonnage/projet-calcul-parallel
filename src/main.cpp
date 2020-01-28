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
static void disp_vector(Vector const& u, int Ncol, int Nlime) {
  cout << "{" << endl;
  for(int i = 0; i < Nlime; ++i) {
    for(int j = 0; j < Ncol; ++j) {
      cout << u[i*Ncol + j] << " ";
    }
    cout << endl;
  }
  cout << "}" << endl;
}

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
        iBegin = iBegin - 1;
    }

    std::cout << "[ " << rank << " ] " << str( iBegin ) " " << iBegin << " " str( iEnd ) " " << iEnd << std::endl;
    //iEnd++;
    int nbLignes = iEnd - iBegin + 1;
    int nbElts   = nbLignes * inputData.colCount;

    float dX = inputData.Lrow / (inputData.rowCount+1), dY = inputData.Lcol / (inputData.colCount+1);
    int   imax = inputData.tMax / inputData.dt;

    float  a = ( ( 2. / ( dX * dX ) ) + ( 2. / ( dY * dY ) ) ) * inputData.D * inputData.dt + 1;
    float  b = (-1. / ( dX * dX )) * inputData.D * inputData.dt;
    float  c = (-1. / ( dY * dY )) * inputData.D * inputData.dt;
    Sparse A( nbLignes, nbLignes, a, b, c );
    std::cout << a << " : " << b << " : " << c << std::endl;

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
    fsource( rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, 0, termeSource, inputData.mode );

    Vector F( nbElts );
    F.set_value( 0. );
    calcul_second_membre( F,
                          nbLignes,
                          inputData.colCount,
                          dX,
                          dY,
                          inputData.D,
                          inputData.dt,
                          gme,
                          hme,
                          termeSource ); // secondMembre dans second_memebre_sparse.f90

    MPI_Datatype line;
    MPI_Type_contiguous( inputData.colCount, MPI_DOUBLE, &line );
    MPI_Type_commit( &line );

    // boucle principale
    double *torecv = new double[inputData.colCount];

    Vector _b(0);
    Vector zero(nbElts);
    zero.set_value(1.);
    imax = 1;
    for ( int i = 0; i < imax; ++i ) {
      _b = F;
      _b.scale( inputData.dt );
      _b += Uprev;
      double error;
      double error_all;
      int j = 0;
      do {
        GC_sparse_parallel( A, _b, U, Uprev, inputData.kmax, inputData.beta );
        Uprev = U;
        int rnext, rprev;
        if(np == 2)
          break;
        rnext = ( rank + 1 ) % np;
        rprev = ( rank + np - 1 ) % np;

        // First send previous last row
        double *tosend = U.data() + ( nbLignes - 1 ) * inputData.colCount;
        // Everyone send to the next rank and receive from the previous
        MPI_Sendrecv( tosend, 1, line, rnext, 0, torecv, 1, line, rprev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        // Except 0, all update gme
        if ( rank != 0 ) { memcpy( gme.data(), torecv, inputData.colCount * sizeof( double ) ); }
        // Second row
        tosend = U.data();
        MPI_Sendrecv( tosend, 1, line, rprev, 1, torecv, 1, line, rnext, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        if ( rank != np - 1 ) {
          memcpy( gme.data() + inputData.colCount, torecv, inputData.colCount * sizeof( double ) );
        }

        error = max_array(U.data(), gme.data(), inputData.colCount);

        for (int i = 0; i < inputData.colCount; i++) {
          cout << U[i] << " ";
        } cout << endl;

        for (int i = 0; i < inputData.colCount; i++) {
          cout << gme[i] << " ";
        } cout << endl;

        MPI_Allreduce(&error, &error_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        cout << rank << " : " << error_all << " : " << i << endl;
        if(rank != 0) {
          memcpy(U.data(), gme.data(), inputData.colCount*sizeof(double));
        }
        if (rank != np - 1) {
          memcpy(U.data() + ( nbLignes - 1 ) * inputData.colCount,
                 gme.data() + inputData.colCount,
                 inputData.colCount*sizeof(double));
        }
        j++;
      } while ( error_all >= inputData.eps && j < 10);
    }
    /*Vector sol(nbElts);
    solution(rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, 0, sol, inputData.mode);

    Vector checker(nbElts);
    checker = sol;
    checker -= U;
    double diff_rel = checker.nrm2();// sol.nrm2();
    std::cout << " diff : " << diff_rel << std::endl;*/

    delete[] torecv;
    MPI_Type_free( &line );
    cout << gme << endl;
    cout << hme << endl;

    //disp_vector(U, inputData.colCount, nbLignes);

    write_vector_to_file( U, inputData.colCount, iBegin, iEnd, dX, dY );
    MPI_Finalize();
    return 0;
}
