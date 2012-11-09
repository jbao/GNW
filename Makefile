# $Id: Makefile 193 2010-10-19 09:35:47Z jbao $

CC                 = g++
CFLAGS             = -Wall -g -fopenmp 
INCLUDES = -I/home/jbao/tool/libsbml-4.2.0/include/ -I/home/jbao/tool/getopt_pp/
LIBDIR     = -L/home/jbao/tool/libsbml-4.2.0/lib/ -L/home/jbao/tool/getopt_pp/
#INSTALL_LIB        = nr2c

#INCLUDEDIR = $(INSTALL_INCLUDEDIR)
#INCLUDES   = -I$(INSTALL_INCLUDEDIR)
#LIBDIR     = -L$(INSTALL_LIBDIR)

LIBS       = -lsbml -lm -lgsl -lgslcblas -lgetopt_pp

#### What should I write here ????????? #### 

TARGET = gnw

OBJS    = GeneNetwork.o HillGene.o Perturbation.o\
	  RandomParameterGaussian.o RandomParameterUniform.o \
	  GnwSettings.o \
	  BenchmarkGenerator.o Edge.o Experiment.o \
	  SteadyStateExperiment.o TimeSeriesExperiment.o \
	  PerturbationMultifactorial.o main.o

$(TARGET): $(OBJS)
		$(CC) $(CFLAGS) $(OBJS) $(LIBDIR) $(LIBS) $(INCLUDES) -o $(TARGET)

.cpp.o:
		$(CC) $(CFLAGS) $(INCLUDES) -c $< 

clean:
		rm -f $(OBJS) $(TARGET)
