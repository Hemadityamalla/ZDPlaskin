FC:=gfortran
RM:=rm -f
OBJFILES:= dvode_f90_m.o example1.o

PROG:= example1
OUTPUT:= *.dat

.PHONY: all clean run


all: $(PROG)
clean: 
	$(RM) $(OBJFILES) *.mod $(PROG) $(OUTPUT) *.log

%.o: %.F90
	gfortran -c $^

$(PROG): $(OBJFILES)
	gfortran -o $@ $^
