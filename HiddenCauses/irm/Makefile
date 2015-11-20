#CC= gnucnoansi
CC= gcc
CC2= gcc

# if compiling with the GNUSL on on MIT's Athena system
# LIBFLAGS = -I/mit/gnusl_v1.5/include
# LINKFLAGS = -L/mit/gnusl_v1.5/lib -lgsl -lgslcblas -lm

LINKFLAGS = -lm

# for debugging with gdb
#CFLAGS= -O0  -g $(LIBFLAGS) 

CFLAGS= -O3 $(LIBFLAGS) 
FINALFLAGS= $(CFLAGS) $(LINKFLAGS)
LDFLAGS= 

LOBJECTS= irm.o irmutils.o chain.o chaincolln.o class.o cokus.o domain.o frontlist.o  item.o multidarray.o opt.o parameters.o relation.o 

LSOURCE=  irm.c irmutils.c chain.c chaincolln.c class.c cokus.c domain.c frontlist.c item.c multidarray.c opt.c parameters.c relation.c 

LHEADERS =  irmutils.h chain.h chaincolln.h class.h cokus.h config.h domain.h frontlist.h item.h multidarray.h opt.h parameters.h relation.h types.h

EXAMPLES = smallbp.opt smallbpmiss.opt smallbpcont.opt smallbpfreqmiss.opt smallbpfreq.opt ring.opt uml.opt  cont1.opt data output/ test/ runtest 

OSHDIR=/afs/athena.mit.edu/user/c/k/ckemp/Experiments/featurelist/results/
OSHFILE=50animalbindat.mat

KINDIR=/afs/csail.mit.edu/u/c/ckemp/RawData/Kinship/Alyawarra/matlabprocess/
KINFILE=alyawarradata.mat

DNATIONSDIR=/afs/csail.mit.edu/u/c/ckemp/RawData/dimensionalitynations/
DNATIONSFILE=dnations.mat


irm:	$(LOBJECTS)
	$(CC) $(FINALFLAGS) $(LOBJECTS) -lm -o irm $(LDFLAGS) -lm

# cokus.c violates some of the conditions gnuc looks for
cokus.o: cokus.c
	$(CC) $(CFLAGS) -c -o cokus.o cokus.c 

# so does opt.c
opt.o: opt.c
	$(CC2) $(CFLAGS) -c -o opt.o opt.c 
tar: 
	tar cvf irmlocal.tar Makefile $(LOCAL) 

tarexternal:
	tar cvf irm.tar Makefile README $(LSOURCE) $(LHEADERS)  $(EXAMPLES) 
	gzip irm.tar

dataexternal:
	tar cvf irmdata.tar README_DATASETS -C data uml.mat \
	-C $(OSHDIR) $(OSHFILE) \
	-C $(KINDIR) $(KINFILE) \
	-C $(DNATIONSDIR) $(DNATIONSFILE) 
	gzip irmdata.tar

clean:
	-rm -f *.o

test:
	./runtest
