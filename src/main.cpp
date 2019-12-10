#include "fonctions"
#include "gradient_conj.h"
#include "second_membre.h"
#include "Input.hpp"
#include "functions.h"

#include <mpi.h>

int main()
{


    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //read_input( rowCount, inputData.Ly, (float)lX, (float)lY /*taille r�elle */, (float)coeff, (float)dt, (float)tmax );
    auto inputData = ReadInput("parametres.dat");

        // param�tres globaux
    // float beta; // pr�cision
    // int       kmax;
    float     dX = inputData.Lx / inputData.Nx, dY = inputData.Ly / inputData.Ny;
    int       imax = inputData.tMax / inputData.dt;

    //setup_procs(); // TODO r�partir les donn�es entre les procs, MPI


    int iBegin, iEnd;
    charge( inputData.Lx,
                     np,
                     rank,
                     iBegin,
                     iEnd ); // retourne l'indice de la ligne de d�but et de celle de fin pour chaque proc

    if(rank != 0){   // toujours avoir une ligne en commun pour tester schwarz
      iBegin = iBegin -1;
    }

    int nbLignes = iEnd - iBegin + 1;
    int nbElts  = nbLignes * inputData.Ly;

    float a, b, c;                                  // les trois valeurs de la matrice penta-diagonale du pdf
    Mat   A = FillA( nbLignes, inputData.Ly, a, b, c ); // Sparse !!! CONST CONST

    Vec U( nbElts ), Uprev( nbElts, 0. ); // Uprev = (0,..., 0)

    /* g : conditions aux bords en haut et en bas (les vecteurs concern�s par MPI !!!!!),
       h : """""""""""""""""""" � gauche et � droite */
    Mat g = init_g( 2, inputData.Ly ); // go fonctions.f90, fonction g
    Mat h = init_h( nbLignes, 2 ); // idem

    Mat termeSource =
        init_terme_source( nbLigne, inputData.Ly ); // CONST CONST. Regarder finstaper dans fonctions.f90 pour l'initialisation

    Vec F( nbELts );
    calcul_second_membre( F, g, h, termeSource ); // secondMembre dans second_memebre_sparse.f90

    /* boucle principale */

    for ( int i = 0; i < imax; ++i ) {
        
        gradient_conj( A, dt * F + Uprev, Uprev, U, beta, kmax, nbElts ); // Go pdf pour comprendre
        Uprev = U;

        communicate_g();// MPI !!!

        calcul_second_membre( F, g, h, termeSource );
    }
}