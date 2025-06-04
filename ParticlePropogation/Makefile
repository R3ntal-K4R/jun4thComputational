
IDIR=include
CC=gcc
CFLAGS=-std=c99 -g -I$(IDIR) 	#Compiler flags used in April 2023 runs
# CFLAGS= -g -I$(IDIR) 			#Compiler flags commonly used in collaborator laptop development

SDIR=src
ODIR=$(SDIR)/obj

LIBS=-lm

_DEPS = auxiliary.h init.h magnetic.h energy_loss.h motion.h generator.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJ = auxiliary.o init.o magnetic.o energy_loss.o motion.o generator.o particleTrajectory.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

particlePropagation: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
