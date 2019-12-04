#include "fonctions"
#include "gradient_conj.h"
#include "second_membre.h"

#include <mpi.h>

int main()
{
    setup_procs(); // répartir les données entre les procs, MPI

    read_input( rowCount, colCount, (float)lX, (float)lY /*taille réelle */, (float)coeff, (float)dt, (float)tmax );

    rank = MPI_COMM_RANK;
    np   = MPI_COMM_SIZE

        // paramètres globaux
        float beta; // précision
    int       kmax;
    float     dX = lX / colCount, dY = lY / rowCount;
    int       imax = tmax / dt;

    int iBegin, iEnd;
    repartir_charge( Nli,
                     np,
                     rank,
                     &iBegin,
                     &iEnd ); // retourne l'indice de la ligne de début et de celle de fin pour chaque proc

    int nbLigne = iEnd - iBegin + 1;
    int nbElts  = nbLignes * colCount;

    float a, b, c;                                  // les trois valeurs de la matrice penta-diagonale du pdf
    Mat   A = FillA( nbLignes, colCount, a, b, c ); // Sparse !!! CONST CONST 

    Vec U( nbElts ), Uprev( nbElts, 0. ); // Uprev = (0,..., 0)

    /* g : conditions aux bords en haut et en bas (les vecteurs concernés par MPI !!!!!),
       h : """""""""""""""""""" à gauche et à droite */
    Mat g = init_g( 2, colCount ); // go fonctions.f90, fonction g
    Mat h = init_h( nbLignes, 2 ); // idem

    Mat termeSource =
        init_terme_source( nbLigne, colCount ); // CONST CONST. Regarder finstaper dans fonctions.f90 pour l'initialisation

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
