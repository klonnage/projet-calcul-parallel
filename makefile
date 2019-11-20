# Nom du compilateur
FC = gfortran

# Directives de compilation
FCFLAGS = -g -Og -fcheck=all -Wall

# Sources (l’ordre est important !)
SRCS = fonctions.f90 syslin.f90  second_membre_sparse.f90  main.f90

# Nom de l’exécutable
PROGRAM = exe

# Ne pas éditer en dessous de cette ligne

OBJ=$(SRCS:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod fort. exe
