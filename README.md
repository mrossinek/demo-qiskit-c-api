# Qiskit C API Demo

This is a very simple demo of Qiskit 2.0's new C API.

## Idea

Implement the Jordan-Wigner fermion-to-qubit mapping.
We do this by parsing an FCIDump file and converting the fermionic terms into
Pauli terms.

## How to use this code

1. Compile the C API for Qiskit @ https://github.com/Cryoris/qiskit-terra/tree/c-api-demo

2. symlink the `qiskit.h` and `libqiskit_cext.so` files into `include/` and `lib/`, respectively

3. symlink your Python library into `lib/libpython.dylib` and the _directory_ of Python headers into `include/python/`

4. run `make run`
