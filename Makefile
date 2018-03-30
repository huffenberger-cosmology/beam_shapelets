OBJ:= beam_shapelets.o  decompose.o  shapelet_cartesian.o
INC:= -I.

all : test $(OBJ)

%.o : %.c
	gcc $(INC) -c -o $@ $< 


test : test.c $(OBJ)
	gcc -Wall $(INC) -o $@ $< $(OBJ) -lm -llapack




