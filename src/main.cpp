#include "GradConj.h"
#include "functions.h"
//#include "second_membre.h"
#include "Input.h"
#include "Sparse.h"
#include "Vector.h"
#include "functions.h"

#include <mpi.h>

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

    if ( rank != 0 ) { // toujours avoir une ligne en commun pour tester schwarz
        iBegin = iBegin - 1;
    }

    int nbLignes = iEnd - iBegin + 1;
    int nbElts   = nbLignes * inputData.Lcol;

    //float a = ((2./(dX*dX))+(2./(dY*dY))*inputData.D*inputData.dt + 1;
    //float b = -1./(dX*dX)*inputData.D*inputData.dt;
    //float c = -1./(dY*dY)*inputData.D*inputData.dt;
    //Sparse A(nbLigne, inputData.Ny, a, b, c);

    //Vector U( nbElts ), Uprev( nbElts );
    //Uprev.set_value(0.);

    ///* g : conditions aux bords en haut et en bas (les vecteurs concern�s par MPI !!!!!),
    //   h : """""""""""""""""""" � gauche et � droite */
    //Vector gme(2 * inputData.Ly);
    //Vector hme(nbLignes * 2);
    //Vector termeSource(nbLignes * inputData.Ly);

    //g(rank, inputData.Ny, dX, inputData.Ly, gme, 1);
    //h(rank, nbLignes, iBegin, dY, inputData.Lx, hme, 1);
    //fsource(rank, inputData.Ny, nbLignes, iBegin, dX, dY, inputData.Lx, inputData.Ly, 0, termeSource, 1);
    //
    //Vector F( nbELts );
    //calcul_second_membre( F, g, h, termeSource ); // secondMembre dans second_memebre_sparse.f90

    ///* boucle principale */

    //for ( int i = 0; i < imax; ++i ) {
    //    gradient_conj( A, dt * F + Uprev, Uprev, U, beta, kmax, nbElts ); // Go pdf pour comprendre
    //    Uprev = U;

    //    communicate_g(); // MPI !!!

    //    calcul_second_membre( F, g, h, termeSource );
    //}
}
