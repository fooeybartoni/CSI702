
EXE1=explicit_central_difference.x
OBJ1=explicit_central_difference.o 

CMP=gcc -g -pg -c -fopenmp
LNK=gcc -pg -fopenmp

all: $(EXE1) 
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE1): $(OBJ1) 
	$(LNK) -o $@ $(OBJ1)

pic.o: explicit_central_difference.c
	$(CMP) $?

