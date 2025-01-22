# Qiskit C API Demo

This is a very simple demo of Qiskit 2.0's new C API.

## Idea

Implement the Jordan-Wigner fermion-to-qubit mapping.
We do this by parsing an FCIDump file and converting the fermionic terms into
Pauli terms.

## How to use this code

1. Compile the C API for Qiskit (obtained from Julien Gacon's branch for the time being)

2. symlink the `qiskit.h` and `libqiskit_cext.so` files into `include/` and `lib/`, respectively

3. run `make run`
