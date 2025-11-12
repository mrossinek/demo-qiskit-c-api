QISKIT_INC=~/Qiskit/qiskit/dist/c/include
QISKIT_LIB=~/Qiskit/qiskit/dist/c/lib
PYTHON_INC=/opt/homebrew/opt/python@3.11/Frameworks/Python.framework/Versions/3.11/include/python3.11
PYTHON_LIB=/opt/homebrew/opt/python@3.11/Frameworks/Python.framework/Versions/3.11/lib
ACCELERATE=~/Qiskit/qiskit/qiskit/_accelerate.cpython-311-darwin.so

libqiskit.so:
	ln -s ${ACCELERATE} libqiskit.so

build: libqiskit.so
	gcc scripts.c -o cextension.so -shared -fpic \
		-I ${QISKIT_INC} -L . -lqiskit \
		-I ${PYTHON_INC} -L ${PYTHON_LIB} -lpython3.11

run_chem: build
	LD_LIBRARY_PATH=./lib python molecule.py

run_ising: build
	LD_LIBRARY_PATH=./lib python ising.py

clean:
	rm cextension.so
