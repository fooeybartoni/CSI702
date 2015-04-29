
EXE3=hello_threads.x
OBJ3=hello_threads.o 

CMP=gcc -g -pg -c -fopenmp
LNK=gcc -pg -fopenmp

all: $(EXE3)
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE3): $(OBJ3) 
	$(LNK) -o $@ $(OBJ3)

hello_threads.o: hello_threads.c
	$(CMP) $?

