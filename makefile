FC = gfortran
DIR = modules

all: lab4 clean

lab4: lab4.f90 file_reader.o keyboard_reader.o parameters.o
	$(FC) lab4.f90 -o lab4 file_reader.o keyboard_reader.o parameters.o


file_reader.o: $(DIR)/file_reader.f90
	$(FC) -c $(DIR)/file_reader.f90

keyboard_reader.o:
	$(FC) -c $(DIR)/keyboard_reader.f90

parameters.o:
	$(FC) -c $(DIR)/parameters.f90

clean:
	del *.o *.mod