
EXE=gradient.x
OBJ=gradient.o 

CMP=gcc -g -pg -c
LNK=gcc -pg -lm

all: $(EXE)
	@echo "done"

clean:
	rm -f *.o *.x 

$(EXE): $(OBJ) 
	$(LNK) -o $@ $(OBJ)

gradient.o: gradient.c
	$(CMP) $?

