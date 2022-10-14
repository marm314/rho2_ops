# Makefile for RHO2_OPS
#
#CPP = g++
CPP = mpicxx -DHAVE_MPI
H5CPP = h5c++
CPPFLAGS = -O3 
PARALLEL = -fopenmp 
Cln = /bin/rm -rf
###########################################
###########################################
SRC=Input_commands.cpp main.cpp Mathematical_Functions.cpp String_ops.cpp gitver.cpp gauss_quad.cpp sphere_lebedev_rule.cpp legendre_quadrature.cpp 
OBJECTS=Input_commands.o main.o Mathematical_Functions.o String_ops.o gitver.o gauss_quad.o sphere_lebedev_rule.o legendre_quadrature.o

all:
	./gitversion.sh 
	make RHO2_OPS
	make chimpanC
	make psi4intface
	make diracintface
	make tar
 
RHO2_OPS: $(OBJECTS) $(SRC) Makefile README 
	$(CPP) $(CPPFLAGS) $(PARALLEL) $(OBJECTS)  -o RHO2_OPS

%.o: %.cpp   
	$(CPP) $(CPPFLAGS) $(PARALLEL) -c $*.cpp 

clean:
	$(Cln) *.o
	$(Cln) *RHO2_OPS
	$(Cln) *chimpanC
	$(Cln) *psi4intface
	$(Cln) *diracintface
	$(Cln) *diracintface_h5
	$(Cln) *~
	$(Cln) gitver.o gitver.cpp gitver.h
	$(Cln) RHO2_OPS.tar.gz 


chimpanC: String_ops.cpp Input_commands.cpp chimpanC.cpp
	$(CPP) $(CPPFLAGS) String_ops.cpp Input_commands_chimpanC.cpp chimpanC.cpp -o chimpanC

psi4intface: psi4intface.cpp
	$(CPP) $(CPPFLAGS) psi4intface.cpp -o psi4intface

diracintface:  String_ops.cpp Input_commands_diracintface.cpp diracintface.cpp diracintface_h5.cpp 
	$(CPP) $(CPPFLAGS) String_ops.cpp Input_commands_diracintface.cpp diracintface.cpp -o diracintface
	$(H5CPP) $(CPPFLAGS) diracintface_h5.cpp -o diracintface_h5


tar:
	tar -pczf RHO2_OPS.tar.gz *.cpp *.h Makefile README test dm2_hf test_psi4 
