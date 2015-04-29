
EXE=gradient.x
OBJ=gradient.o 

CMP=gcc -g -pg -c -fopenmp
LNK=gcc -pg -fopenmp

all: $(EXE)
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE): $(OBJ) 
	$(LNK) -o $@ $(OBJ)

gradient.o: gradient.c
	$(CMP) $?

