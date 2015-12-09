#
#
F90 = gfortran
LIB = 
DEBUG_FLAG = -g -fbounds-check -fdefault-real-8
OPTIM_FLAG = -O3 -fdefault-real-8
USUAL_FLAG = -fdefault-real-8
PROFIL_FLAG = -g -pg -fdefault-real-8

PROG = run
SRC = display.f90 transport.f90 main.f90

usual :
	$(F90) $(USUAL_FLAG) $(SRC) -o $(PROG) $(LIB)

debug :
	$(F90) $(DEBUG_FLAG) $(SRC) -o $(PROG) $(LIB)

profil:
	$(F90) $(PROFIL_FLAG) $(SRC) -o $(PROG) $(LIB)

optim :
	$(F90) $(OPTIM_FLAG) $(SRC) -o $(PROG) $(LIB)


clean :
	@rm -f *.o *.mod *~ core a.out 
	@echo "On a fait du nettoyage"
