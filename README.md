# Projet de Calcul Parallel

## Compilation

Créer une répertoire build, se placer dedans et lancer la commande cmake puis make:

```sh
mkdir build
cd build
cmake3 ..
make
```

## Exécution

L'exécutable se situe dans le répertoire **src** ainsi créer.
Avant de lancer le programme, il faut se placer dans **src** et créer le répertoire **out**:

```sh
mkdir -p out
```

Pour lancer le programme, se placer dans **src/** et lancer:

```sh
mpirun -np 1 exe
```

Cela lance le programme en séquentiel et écrit la solution dans **out**.

Changer 1 en le nombre de processus souhaiter pour avoir une exécution séquentielle.

## Paramètres du problème

Les paramètres sont dans le fichier **parametres.dat** dans le répertoire **data/**
à la racine.

## Affichage des données

La solution est écrite dans le répertoire **out** placé à l'endroit où est créer l'exécutable.

Pour afficher la solution, il faut lancer `gnuplot` et afficher successivement les courbes.

Le modèle:
```sh
splot "out/GradConj_[rank du processus]_[Jour de la semaine]_[heure].txt" u 1:2:3
```

Un exemple:
```sh
splot "out/GradConj_0_Fri_Jan_31_11_38_10_2020_.txt" u 1:2:3
replot "out/GradConj_1_Fri_Jan_31_11_38_10_2020_.txt" u 1:2:3
```