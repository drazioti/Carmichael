#CC is the compiler used.
CC =  g++
#LIBS used from the compiler.
LIBS =  -lgmpxx -lgmp -lcrypto
LIBS_PAR =  -fopenmp
ALL_LIBS =
#CFLAGS passed to the compiler.
#-fomit-frame-pointer flag does not help much
CFLAGS =  -std=c++11 -c -O3
#Names of main objects
OBJ_MAIN = carmi.o
# Names of all objects exept main
OBJS_NO_MAIN = subset_product.o Combinations.o
# Names of executables
EXECUTABLE_MAIN = carmichael.out
EXECUTABLE_MAIN_PAR = carmichael_par.out
# Names for all files of certain type
CPP_FILES     = $(wildcard *.cpp)
OBJ_FILES     = $(wildcard *.o)
#clean does not work with wildcard when .o files are not compiled
#thats why I use .o on clean_obj

# Makefile commands start here
all: clean
	make serial clean_obj
	make parallel clean_obj
	@echo "Build Carmichael executable"
	
serial: ALL_LIBS := $(LIBS)
serial: $(OBJ_MAIN) $(OBJS_NO_MAIN)
	$(CC) $(OBJ_MAIN) $(OBJS_NO_MAIN) -o $(EXECUTABLE_MAIN) $(ALL_LIBS)

parallel: ALL_LIBS := $(LIBS) $(LIBS_PAR)
parallel: $(OBJ_MAIN) $(OBJS_NO_MAIN)
	$(CC) $(OBJ_MAIN) $(OBJS_NO_MAIN) -o $(EXECUTABLE_MAIN_PAR)  $(ALL_LIBS)

prepare: $(OBJ_MAIN) $(OBJS_NO_MAIN)

carmi.o:  carmi.cpp
	$(CC) $(CFLAGS) carmi.cpp $(ALL_LIBS)

subset_product.o:  subset_product.cpp
	$(CC) $(CFLAGS) subset_product.cpp  $(ALL_LIBS)

Combinations.o:  Combinations.cpp
	$(CC) $(CFLAGS) Combinations.cpp $(ALL_LIBS)

clean: clean_obj
	rm -f $(EXECUTABLE_MAIN)
	rm -f $(EXECUTABLE_MAIN_PAR)
	@echo "Cleanned executables"

clean_obj:
	rm -f *.o
	@echo "Cleanned objs"
