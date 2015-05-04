### Basic Makefile for MPI

CC      = mpicc
CFLAGS  = -g -O0
LD      = mpicc
LDFLAGS = -g

PROGRAM = hello-mpi

all:	${PROGRAM}.exe

${PROGRAM}.exe:			${PROGRAM}.o
		${LD} ${LDFLAGS} $< -o ${PROGRAM}.exe

${PROGRAM}.o:           ${PROGRAM}.c
		${CC} ${CFLAGS} -c $< -o ${PROGRAM}.o

clean:
		rm -f ${PROGRAM}.o ${PROGRAM}.exe