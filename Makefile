OUTPUT_FOLDER = output

all: run_hexagonal

run_hexagonal:
	$(MAKE) -C mesh/ hexagonal.msh
	python3 demo.py --mesh mesh/hexagonal.msh 1 11.8

run_hexagonal_and_plot:
	$(MAKE) -C mesh/ hexagonal.msh
	python3 demo.py --mesh mesh/hexagonal.msh -o $(OUTPUT_FOLDER) 1 11.8

run_hexagonal_and_log:
	$(MAKE) -C mesh/ hexagonal.msh
	python3 demo.py --mesh mesh/hexagonal.msh 1 11.8 -log_view

clean:
	$(MAKE) -C mesh/ clean
