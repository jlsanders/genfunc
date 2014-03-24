
# Makefile for new_struct

#------------------------------------------------------
CCOMPILER = g++-4.8 -std=c++0x #clang++ -std=c++11 -stdlib=libc++ -nostdinc++ -I/usr/include/c++/4.2.1/
BOOSTINCS = -I/usr/local/Cellar/boost/1.55.0/include/boost/
BOOSTLIBS = -L/usr/local/Cellar/boost/1.55.0/lib/ -lboost_python-mt -lboost_python
PYTHONLIBS = -L/usr/local/Cellar/python/2.7.3/Frameworks/Python.framework/Versions/2.7/lib/ -lpython2.7

OPTS = -Wall -O3 -ffast-math -fno-strict-aliasing
INC = -Iinc/ -I../../ $(BOOSTINCS)
SRC = src/
POT = pot/
OBJ = obj/
LIB	= $(BOOSTLIBS) $(PYTHONLIBS) -lgsl -lgslcblas -lm

ALL = triax_py.so
OBJECTS = $(OBJ)coordsys.o $(OBJ)coordtransforms.o $(OBJ)potential.o $(OBJ)orbit.o $(OBJ)aa.o

#------------------------------------------------------

default: clean all

all:  $(ALL)

forceall:	clean all

#------------------------------------------------------

obj/%.o : src/%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

#------------------------------------------------------

triax_py.so: $(OBJECTS) $(OBJ)triax_py.o 
	$(CCOMPILER) -fPIC -shared $(CFL) $(OBJECTS) $(OBJ)triax_py.o $(INC) $(LIB) $(LIBS) -o triax_py.so

clean:
	rm -rf *~; rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe $(ALL)
