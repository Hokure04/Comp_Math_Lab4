FC = gfortran
DIR = modules
LOG = logarithmic_approximation
POW = power_approximation
EXP = exponential_approximation
CUBE = cubic_approximation
QUADR = quadratic_approximation
LINE = linear_approximation
PARAM = parameters
KEY = keyboard_reader
FILE = file_reader

all: lab4 clean

lab4: lab4.f90 $(FILE).o $(KEY).o $(PARAM).o $(LINE).o $(QUADR).o $(CUBE).o $(EXP).o $(POW).o $(LOG).o solve_system.o display_graphics.o
	$(FC) lab4.f90 -o lab4 $(FILE).o $(KEY).o $(PARAM).o $(LINE).o $(QUADR).o $(CUBE).o $(EXP).o $(POW).o $(LOG).o solve_system.o display_graphics.o

file_reader.o: $(DIR)/$(FILE).f90
	$(FC) -c $(DIR)/$(FILE).f90

keyboard_reader.o: $(DIR)/$(KEY).f90
	$(FC) -c $(DIR)/$(KEY).f90

display_graphics.o: $(DIR)/display_graphics.f90
	$(FC) -c $(DIR)/display_graphics.f90

parameters.o: $(DIR)/$(PARAM).f90
	$(FC) -c $(DIR)/$(PARAM).f90

linear_approximation.o: $(DIR)/$(LINE).f90
	$(FC) -c $(DIR)/$(LINE).f90

quadratic_approximation.o: $(DIR)/$(QUADR).f90
	$(FC) -c $(DIR)/$(QUADR).f90

cubic_approximation.o: $(DIR)/$(CUBE).f90
	$(FC) -c $(DIR)/$(CUBE).f90

solve_system.o: $(DIR)/solve_system.f90
	$(FC) -c $(DIR)/solve_system.f90

exponential_approximation.o: $(DIR)/$(EXP).f90
	$(FC) -c $(DIR)/$(EXP).f90 

power_approximation.o: $(DIR)/$(POW).f90
	$(FC) -c $(DIR)/$(POW).f90

logarithmic_approximation.o: $(DIR)/$(LOG).f90
	$(FC) -c $(DIR)/$(LOG).f90

clean:
	del *.o *.mod