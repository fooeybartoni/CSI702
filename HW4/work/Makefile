
EXE=mpi_hw4a.x
OBJ=mpi_hw4a.o 

CMP=gcc -g -pg -c
LNK=gcc -pg -lm

all: $(EXE)
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE): $(OBJ) 
	$(LNK) -o $@ $(OBJ)

gradient.o: mpi_hw4a.c
	$(CMP) $?

