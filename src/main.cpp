#include "GradConj.h"
#include "Input.h"
#include "Sparse.h"
#include "Vector.h"
#include "functions.h"
#include "second_membre.h"

#include <cstring>
#include <iostream>
#include <mpi.h>

#define str( x ) #x

int main( int argc, char **argv )
{
    int rank, np;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &np );

    auto inputData = ReadInput( "../../data/parametres.dat" );

    float dX = inputData.Lrow / inputData.rowCount, dY = inputData.Lcol / inputData.colCount;
    int   imax = inputData.tMax / inputData.dt;

    // setup_procs(); // TODO r�partir les donn�es entre les procs, MPI

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

    int nbLignes = iEnd - iBegin + 1;
    int nbElts   = nbLignes * inputData.Lcol;

    float  a = ( 2. / ( dX * dX ) ) + ( 2. / ( dY * dY ) ) * inputData.D * inputData.dt + 1;
    float  b = -1. / ( dX * dX ) * inputData.D * inputData.dt;
    float  c = -1. / ( dY * dY ) * inputData.D * inputData.dt;
    Sparse A( nbLignes, inputData.rowCount, a, b, c );

    Vector U( nbElts ), Uprev( nbElts );
    Uprev.set_value( 0. );

    ///* g : conditions aux bords en haut et en bas (les vecteurs concern�s par MPI !!!!!),
    //   h : """""""""""""""""""" � gauche et � droite */
    Vector gme( 2 * inputData.Lcol );
    Vector hme( nbLignes * 2 );
    Vector termeSource( nbLignes * inputData.Lcol );

    g( rank, inputData.colCount, dX, inputData.Lcol, gme, 1 );
    h( rank, nbLignes, iBegin, dY, inputData.Lrow, hme, 1 );
    fsource( rank, inputData.colCount, nbLignes, iBegin, dX, dY, inputData.Lrow, inputData.Lcol, 0, termeSource, 1 );

    Vector F( nbElts );
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

    ///* boucle principale */
    double *torecv = new double[inputData.colCount];

    for ( int i = 0; i < imax; ++i ) {
        //    gradient_conj( A, dt * F + Uprev, Uprev, U, beta, kmax, nbElts ); // Go pdf pour comprendre
        //    Uprev = U;

        Vector b = F;
        b.scale( inputData.dt );
        b += Uprev;
        double error;
        do {
            GC_sparse_parallel( A, b, U, Uprev, inputData.kmax, inputData.beta );
            Uprev = U;
            int rnext, rprev;
            rnext = ( rank + 1 ) % np;
            rprev = ( rank + np - 1 ) % np;

            /* First send previous last row */
            double *tosend = U.data() + ( nbLignes - 2 ) * inputData.colCount;
            /* Everyone send to the next rank and receive from the previous */
            MPI_Sendrecv( tosend, 1, line, rnext, 0, torecv, 1, line, rprev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            /* Except 0, all update gme */
            if ( rank != 0 ) { memcpy( gme.data(), torecv, inputData.colCount * sizeof( double ) ); }
            /* Second row */
            tosend = U.data() + inputData.colCount;
            MPI_Sendrecv( tosend, 1, line, rprev, 0, torecv, 1, line, rnext, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
            if ( rank != np - 1 ) {
                memcpy( gme.data() + inputData.colCount, torecv, inputData.colCount * sizeof( double ) );
            }
        } while ( error >= inputData.eps );
    }

    delete[] torecv;
    MPI_Type_free( &line );
    MPI_Finalize();
    return 0;
}
