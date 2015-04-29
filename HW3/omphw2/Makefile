
EXE3=pic.x
OBJ3=pic.o 

CMP=gcc -g -pg -c -fopenmp
LNK=gcc -pg -fopenmp
FCM=gfortran -c -fopenmp
FLN=gfortran -fopenmp

all: $(EXE1) $(EXE2) $(EXE3)
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE3): $(OBJ3) 
	$(LNK) -o $@ $(OBJ3)

pic.o: pic.c
	$(CMP) $?

