all: 
	gfortran run_model.f90 basic.f90 tools.f90 model_tools.f90 -o model.out
	./model.out 
clean:
	rm *.out