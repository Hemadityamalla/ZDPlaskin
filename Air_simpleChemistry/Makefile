FC:=gfortran
RM:=rm -f
OBJFILES:= dvode_f90_m.o zdplaskin_m.o main.o

LINKLIB:= bolsig_x86_64.so

PROG:= simpleChem
OUTPUT:= qt*.txt

.PHONY: all clean run


all: $(PROG)

run: $(PROG)
	./$(PROG)

clean: 
	$(RM) $(OBJFILES) *.mod $(PROG) $(OUTPUT) *.log

%.o: %.F90
	$(FC) -c $^

$(PROG): $(OBJFILES) $(LINKLIB)
	$(FC) -o $@ $^


zdplaskin_m.F90: kinetics.inp preprocessor
	./preprocessor kinetics.inp $@
	


