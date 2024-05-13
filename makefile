FC = gfortran
DIR = modules


all: lab4 clean

lab4: lab4.f90 file_reader.o keyboard_reader.o parameters.o linear_aproximation.o quadratic_aproximation.o cubic_aproximation.o exponential_approximation.o power_approximation.o
	$(FC) lab4.f90 -o lab4 file_reader.o keyboard_reader.o parameters.o linear_aproximation.o quadratic_aproximation.o cubic_aproximation.o exponential_approximation.o power_approximation.o

file_reader.o: $(DIR)/file_reader.f90
	$(FC) -c $(DIR)/file_reader.f90

keyboard_reader.o: $(DIR)/keyboard_reader.f90
	$(FC) -c $(DIR)/keyboard_reader.f90

parameters.o: $(DIR)/parameters.f90
	$(FC) -c $(DIR)/parameters.f90

linear_aproximation.o: $(DIR)/linear_aproximation.f90
	$(FC) -c $(DIR)/linear_aproximation.f90

quadratic_aproximation.o: $(DIR)/quadratic_aproximation.f90
	$(FC) -c $(DIR)/quadratic_aproximation.f90

cubic_aproximation.o: $(DIR)/cubic_aproximation.f90
	$(FC) -c $(DIR)/cubic_aproximation.f90

exponential_approximation.o: $(DIR)/exponential_approximation.f90
	$(FC) -c $(DIR)/exponential_approximation.f90

power_approximation.o: $(DIR)/power_approximation.f90
	$(FC) -c $(DIR)/power_approximation.f90

clean:
	del *.o *.mod