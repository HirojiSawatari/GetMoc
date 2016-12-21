FC = ifort
FFLAGS = -O2 
INC = -I/fs/local/include
LIBS = -L/fs/local/lib -lnetcdf
TARGETS = pop_moc 

all: $(TARGETS)

pop_moc : pop_moc.f90 
	$(FC) $(FFLAGS)  -o pop_moc pop_moc.f90 $(INC) $(LIBS) 
pop_moc_ei : pop_moc_ei.f90 
	$(FC) $(FFLAGS)  -o pop_moc_ei pop_moc_ei.f90 $(INC) $(LIBS) 
#pop_moc : basins_pop_moc.f90 
#	$(FC) $(FFLAGS)  -o pop_moc basins_pop_moc.f90 $(INC) $(LIBS) 
