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

#define str( x ) #x

using namespace std;
static void disp_vector(Vector const& u, int Ncol) {
  cout << "{" << endl;
  for(int i = 0; i < Ncol; ++i) {
    for(int j = 0; j < Ncol; ++j) {
      cout << u[i*Ncol + j] << " ";
    }
    cout << endl;
  }
  cout << "}" << endl;
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

    std::cout << "[ " << rank << " ] " << str( iBegin ) " " << iBegin << " " str( iEnd ) " " << iEnd << std::endl;

    if ( rank != 0 ) { // toujours avoir une ligne en commun pour tester schwarz
        iBegin = iBegin - 1;
    }

    iEnd++;
    int nbLignes = iEnd - iBegin + 1;
    int nbElts   = nbLignes * inputData.colCount;

    float dX = inputData.Lrow / (inputData.rowCount+1), dY = inputData.Lcol / (inputData.colCount+1);
    int   imax = inputData.tMax / inputData.dt;

    float  a = ( ( 2. / ( dX * dX ) ) + ( 2. / ( dY * dY ) ) ) * inputData.D * inputData.dt + 1;
    float  b = (-1. / ( dX * dX )) * inputData.D * inputData.dt;
    float  c = (-1. / ( dY * dY )) * inputData.D * inputData.dt;
    Sparse A( inputData.colCount, inputData.colCount, a, b, c );
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
    for ( int i = 0; i < imax; ++i ) {
      _b = F;
      _b.scale( inputData.dt );
      _b += Uprev;
      double error = 0.;
      do {
        GC_sparse_parallel( A, _b, U, Uprev, inputData.kmax, inputData.beta );
        
        Uprev = U;
        int rnext, rprev;
        if(np == 1)
          break;
        rnext = ( rank + 1 ) % np;
        rprev = ( rank + np - 1 ) % np;

        // First send previous last row
        double *tosend = U.data() + ( nbLignes - 2 ) * inputData.colCount;
        // Everyone send to the next rank and receive from the previous
        MPI_Sendrecv( tosend, 1, line, rnext, 0, torecv, 1, line, rprev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        // Except 0, all update gme
        if ( rank != 0 ) { memcpy( gme.data(), torecv, inputData.colCount * sizeof( double ) ); }
        // Second row
        tosend = U.data() + inputData.colCount;
        MPI_Sendrecv( tosend, 1, line, rprev, 0, torecv, 1, line, rnext, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        if ( rank != np - 1 ) {
          memcpy( gme.data() + inputData.colCount, torecv, inputData.colCount * sizeof( double ) );
        }
      } while ( error >= inputData.eps );
    }
    Vector sol(nbElts);
    solution(rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, 0, sol, inputData.mode);

    Vector checker(nbElts);
    checker = sol;
    checker -= U;
    double diff_rel = checker.nrm2();// sol.nrm2();
    std::cout << " diff : " << diff_rel << std::endl;

    /*std::cout << "============== Display vector U ==============" << std::endl;

    disp_vector (U, inputData.rowCount);

    F.set_value(1.);
    F.scale(inputData.dt);
    U = A*F;

    std::cout << "============== Display vector sol ==============" << std::endl;
    std::cout << a << " : " << b << " : " << c << std::endl;  
    /*Sparse B(inputData.rowCount, inputData.colCount, a, b, c);
    GC_sparse_parallel(B, F, U, zero, inputData.kmax, inputData.beta);    
    disp_vector (sol, inputData.rowCount);

    std::cout << "============== My checker ==============" << std::endl;
    GC_sparse_parallel(A, F, U, zero, inputData.kmax, inputData.beta/10);    
    Vector prod;
    prod = A*U;
    prod -= F;
    double diff_sol = prod.nrm2();
    std::cout << " real sol : " << diff_sol << std::endl;*/

    delete[] torecv;
    MPI_Type_free( &line );

    write_vector_to_file( U, inputData.colCount, iBegin, iEnd, dX, dY );
    MPI_Finalize();
    return 0;
}
