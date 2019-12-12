FC:=gfortran
FFLAGS:= -w
RM:=rm -f
OBJFILES:= dvode_f90_m.o zdplaskin_m.o main.o

LINKLIB:= bolsig_x86_64.so

PROG:= AirTest
OUTPUT:= *.txt

.PHONY: all clean run


all: $(PROG)

run: $(PROG)
	./$(PROG)

clean: 
	$(RM) $(OBJFILES) *.mod *.log $(PROG) $(OUTPUT)

%.o: %.F90
	$(FC) $(FFLAGS) -c $^

$(PROG): $(OBJFILES) $(LINKLIB)
	$(FC) -o $@ $^


#We do not use the below commands as we do not (yet) have the .inp file for air mechanisms-- work for later

#zdplaskin_m.F90: kinetics.inp preprocessor
#	./preprocessor kinetics.inp $@
	


