// for Python interface
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// main C script
#include "qiskit.h"
#include <assert.h>
#include <complex.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

QkSparseObservable *jw_term(uintmax_t num_qubits, uintmax_t site, complex double coeff, bool creation)
{
  QkSparseObservable *op = qk_obs_zero(num_qubits);

  uint32_t indices[site];
  QkBitTerm x_bits[site];
  QkBitTerm y_bits[site]; // could just allocate 1 of these

  for (uintmax_t i = 0; i < site - 1; i++)
  {
    indices[i] = i;
    x_bits[i] = QkBitTerm_Z;
    y_bits[i] = QkBitTerm_Z;
  }
  indices[site - 1] = site - 1;
  x_bits[site - 1] = QkBitTerm_X;
  y_bits[site - 1] = QkBitTerm_Y;

  complex double x_coeff = 0.5 * coeff;
  complex double y_coeff = coeff * (creation ? -0.5 * I : 0.5 * I);

  QkSparseTerm x = {x_coeff, site, x_bits, indices, num_qubits};
  QkSparseTerm y = {y_coeff, site, y_bits, indices, num_qubits};

  qk_obs_add_term(op, &x);
  qk_obs_add_term(op, &y);

  return op;
}

QkSparseObservable *get_qubit_observable()
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen("h2.fcidump", "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  uintmax_t num_orbs;
  uintmax_t num_qubits;

  QkSparseObservable *obs;

  int line_co = 0;
  while ((read = getline(&line, &len, fp)) != -1)
  {
    line_co++;
    if (line_co == 1)
    {
      // get the number of orbitals
      sscanf(line, "&FCI NORB = %ju,", &num_orbs);
      num_qubits = 2 * num_orbs;
      obs = qk_obs_zero(num_qubits);
    }

    if (line_co < 5)
    {
      continue;
    }
    char *end = NULL;
    double c = strtod(line, &end);
    double complex coeff = c + 0.0 * I;
    uintmax_t i = strtoumax(end, &end, 10);
    uintmax_t a = strtoumax(end, &end, 10);
    uintmax_t j = strtoumax(end, &end, 10);
    uintmax_t b = strtoumax(end, &end, 10);

    if (i == 0 && j == 0 && a == 0 && b == 0)
    {
      uint32_t inds[] = {};
      QkBitTerm bits[] = {};
      QkSparseTerm energy = {coeff, 0, bits, inds, num_qubits};
      qk_obs_add_term(obs, &energy);
    }
    else if (j == 0 && b == 0)
    {
      for (int is_beta = 0; is_beta <= 1; is_beta++)
      {
        QkSparseObservable *a_i = jw_term(num_qubits, i + is_beta * num_orbs, coeff, true);
        QkSparseObservable *a_a = jw_term(num_qubits, a + is_beta * num_orbs, coeff, false);
        QkSparseObservable *out = qk_obs_compose(a_a, a_i);

        qk_obs_append(obs, out);

        qk_obs_free(a_i);
        qk_obs_free(a_a);
        qk_obs_free(out);
      }
    }
    else
    {
      for (int is_beta = 0; is_beta <= 1; is_beta++)
      {
        // I don't really know what's the correct combination here
        QkSparseObservable *a_i = jw_term(num_qubits, i + is_beta * num_orbs, coeff, true);
        QkSparseObservable *a_j = jw_term(num_qubits, j + is_beta * num_orbs, coeff, true);
        QkSparseObservable *a_a = jw_term(num_qubits, a + is_beta * num_orbs, coeff, false);
        QkSparseObservable *a_b = jw_term(num_qubits, b + is_beta * num_orbs, coeff, false);
        QkSparseObservable *out = qk_obs_compose(a_b, a_a);
        out = qk_obs_compose(out, a_j);
        out = qk_obs_compose(out, a_i);

        qk_obs_append(obs, out);

        qk_obs_free(a_i);
        qk_obs_free(a_j);
        qk_obs_free(a_a);
        qk_obs_free(a_b);
        qk_obs_free(out);
      }
    }
  }

  fclose(fp);
  if (line)
    free(line);

  obs = qk_obs_canonicalize(obs, 1e-10);

  return obs;
}

// build the PyCapsule containing the sparse observable
static PyObject *cmod_qubit_observable(PyObject *self, PyObject *args)
{
  QkSparseObservable *obs = get_qubit_observable();
  PyObject *capsule;
  capsule = PyCapsule_New((void *)obs, "cbuilder.qubit_observable", NULL);
  return capsule;
}

static PyMethodDef CModMethods[] = {
    {"get_qubit_observable", cmod_qubit_observable, METH_VARARGS, "Get the qubit observable"},
    {NULL, NULL, 0, NULL}, // sentinel
};

static struct PyModuleDef cmod = {
    PyModuleDef_HEAD_INIT,
    "cmod", // module name
    NULL,   // docs
    -1,     // keep the module state in global variables
    CModMethods,
};

PyMODINIT_FUNC PyInit_cmod(void) { return PyModule_Create(&cmod); }

int main(int argc, char *argv[])
{
  PyStatus status;
  PyConfig config;
  PyConfig_InitPythonConfig(&config);

  /* Add a built-in module, before Py_Initialize */
  if (PyImport_AppendInittab("cmod", PyInit_cmod) == -1)
  {
    fprintf(stderr, "Error: could not extend in-built modules table\n");
    exit(1);
  }

  /* Pass argv[0] to the Python interpreter */
  status = PyConfig_SetBytesString(&config, &config.program_name, argv[0]);
  if (PyStatus_Exception(status))
  {
    goto exception;
  }

  /* Initialize the Python interpreter.  Required.
     If this step fails, it will be a fatal error. */
  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status))
  {
    goto exception;
  }
  PyConfig_Clear(&config);

  /* Optionally import the module; alternatively,
     import can be deferred until the embedded script
     imports it. */
  PyObject *pmodule = PyImport_ImportModule("cmod");
  if (!pmodule)
  {
    PyErr_Print();
    fprintf(stderr, "Error: could not import module 'cmod'\n");
  }

  return 0;

exception:
  PyConfig_Clear(&config);
  Py_ExitStatusException(status);
}
