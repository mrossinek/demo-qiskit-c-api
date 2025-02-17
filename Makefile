which-python:
	python -c 'from distutils.sysconfig import get_config_var; print("lib:", get_config_var("LIBDIR"))'
	python -c 'from sysconfig import get_path; print("include:", get_path("include"))'

build:
	gcc main.c -o cmod.so -shared -fpic -I./include -I./include/python -lqiskit_cext -lpython -L./lib 

run: build
	LD_LIBRARY_PATH=./lib python main.py

run_compose:
	gcc load_fci.c -o cmod.so -shared -fpic -I./include -I./include/python -lqiskit_cext -lpython -L./lib 
	LD_LIBRARY_PATH=./lib python main.py

run_ising:
	gcc load_fci.c -o cmod.so -shared -fpic -I./include -I./include/python -lqiskit_cext -lpython -L./lib 
	LD_LIBRARY_PATH=./lib python ising.py
