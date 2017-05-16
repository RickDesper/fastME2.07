#revised 12/11/2002
#revised 12/14/2004
#revised 9/6/2005

CC = gcc
EXEC = fastme
LIBS = -lm
#CFLAGS= -g
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -Wall -static
DFLAG = -DUNIX

######################################################################################################

#MEOBJS = graph.o inputs.o gme.o bme.o  distance.o newick.o heap.o traverse.o NNI.o bNNI.o random.o SPR.o TBR.o print.o nj.o FastDist.o  CompDist.o fastme.o interface_utilities.o interface_free.o interface_options.o

MEOBJS = graph.o inputs.o gme.o bme.o  distance.o newick.o heap.o traverse.o NNI.o bNNI.o random.o SPR.o TBR.o FastDist.o  CompDist.o fastme.o interface_utilities.o interface_options.o BIONJ.o utils.o MVR.o p_utils.o p_eigen.o p_models.o p_lk.o p_optimiz.o p_bootstrap.o

 #ambiguity.o

$(EXEC) : fastme.c $(MEOBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(MEOBJS) $(LIBS)

clean:
	@rm *.o
######################################################################################################


#the next command needs only to be called once.  If the directory already 
#exists, delete the command.
#	mkdir ./bin
#BINDIR=./bin
#alternatively, change next line to local binary directory as desired
#BINDIR= ~rick/bin
#fastme: fastme.c $(MEOBJS)
#	$(CC) $(CFLAGS) -o $(BINDIR)/fastme fastme.c $(MEOBJS) $(LIBS) 
fastme.o: fastme.c
	$(CC) $(CFLAGS) $(DFLAG) -c fastme.c
distance.o: distance.c
	$(CC) $(CFLAGS) $(DFLAG) -c distance.c
graph.o: graph.c
	$(CC) $(CFLAGS) $(DFLAG) -c graph.c
inputs.o: inputs.c
	$(CC) $(CFLAGS) $(DFLAG) -c inputs.c
gme.o: gme.c
	$(CC) $(CFLAGS) $(DFLAG) -c gme.c
bme.o: bme.c
	$(CC) $(CFLAGS) $(DFLAG) -c bme.c
nj.o: nj.c
	$(CC) $(CFLAGS) $(DFLAG) -c nj.c
newick.o: newick.c
	$(CC) $(CFLAGS) $(DFLAG) -c newick.c
heap.o: heap.c
	$(CC) $(CFLAGS) $(DFLAG) -c heap.c
traverse.o: traverse.c
	$(CC) $(CFLAGS) $(DFLAG) -c traverse.c
NNI.o: NNI.c
	$(CC) $(CFLAGS) $(DFLAG) -c NNI.c
bNNI.o: bNNI.c
	$(CC) $(CFLAGS) $(DFLAG) -c bNNI.c -o bNNI.o
SPR.o: SPR.c
	$(CC) $(CFLAGS) $(DFLAG) -c SPR.c -o SPR.o
TBR.o: TBR.c
	$(CC) $(CFLAGS) $(DFLAG) -c TBR.c -o TBR.o
#print.o: print.c
#	$(CC) $(CFLAGS) $(DFLAG) -c print.c -o print.o
random.o: random.c
	$(CC) $(CFLAGS) $(DFLAG) -c random.c -o random.o
FastDist.o: FastDist.c
	$(CC) $(CFLAGS) $(DFLAG) -c FastDist.c -o FastDist.o
CompDist.o: CompDist.c
	$(CC) $(CFLAGS) $(DFLAG) -c CompDist.c -o CompDist.o 
interface_utilities.o : interface_utilities.c interface_utilities.h
	$(CC) $(CFLAGS) $(DFLAG) -c interface_utilities.c
#interface_free.o : interface_free.c interface_free.h
#	$(CC) $(CFLAGS) $(DFLAG) -c interface_free.c
interface_options.o : interface_options.c interface_options.h
	$(CC) $(CFLAGS) $(DFLAG) -c interface_options.c
BIONJ.o : BIONJ.c BIONJ.h
	$(CC) $(CFLAGS) $(DFLAG) -c BIONJ.c
utils.o : utils.c utils.h
	$(CC) $(CFLAGS) $(DFLAG) -c utils.c
MVR.o : MVR.c MVR.h
	$(CC) $(CFLAGS) $(DFLAG) -c MVR.c

p_utils.o : p_utils.c p_utils.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_utils.c
p_eigen.o : p_eigen.c p_eigen.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_eigen.c
p_models.o : p_models.c p_models.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_models.c
p_lk.o : p_lk.c p_lk.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_lk.c
p_optimiz.o : p_optimiz.c p_optimiz.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_optimiz.c
p_bootstrap.o : p_bootstrap.c p_bootstrap.h
	$(CC) $(CFLAGS) $(DFLAG) -c p_bootstrap.c


