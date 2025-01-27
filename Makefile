which-python:
	python -c 'from distutils.sysconfig import get_config_var; print("lib:", get_config_var("LIBDIR"))'
	python -c 'from sysconfig import get_path; print("include:", get_path("include"))'

build:
	gcc main.c -o cmod.so -shared -fpic -I./include -I./include/python -lqiskit_cext -lpython -L./lib 

run: build
	python main.py
