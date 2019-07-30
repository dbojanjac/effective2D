#OUTPUT_FOLDER = ...
#PLOT_FOLDER = ...

all: run_hexagonal

run_hexagonal:
	$(MAKE) -C mesh/ hexagonal.msh
	python3 demo.py --mesh mesh/hexagonal.msh 1 11.8

# run python script, get effective params and save pictures from matplotlib
#run_experiment_and_plot:

# get effective params
#effective:

# benchmark running effective with PETSc -log_view
#benchmark:

clean:
	$(MAKE) -C mesh/ clean
