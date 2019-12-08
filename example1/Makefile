FC:=gfortran
RM:=rm -f
OBJFILES:=*.o

LINKLIB:=bolsig_x86_64.so

PROG:=example1
OUTPUT:=out.txt

.PHONY: all clean run

all: $(PROG)

clean: 
	$(RM) $(OBJFILES) *.mod $(PROG) $(OUTPUT)

%.o: %.F90
	$(FC) -c $@ $^

$(PROG): 
