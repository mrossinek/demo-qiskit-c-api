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

    if (i == 0)
    {
      // constant energy offset
      assert(a == 0);
      assert(j == 0);
      assert(b == 0);
      // coeff * Id
      QkBitTerm bit_terms[] = {};
      uint32_t indices[] = {};
      QkSparseTerm term = {coeff, 0, bit_terms, indices, num_qubits};
      qk_obs_add_term(obs, &term);
    }
    else if (j == 0)
    {
      // 1-body term
      assert(b == 0);
      assert(i != 0);
      assert(a != 0);

      if (i == a)
      {
        QkSparseObservable *terms = qk_obs_identity(num_qubits);
        // coeff * Id
        terms = qk_obs_multiply(terms, &coeff);
        // alpha-spin: -0.5 * coeff * Z_i
        double complex coeff_half = -0.5 * coeff;
        QkBitTerm bit_terms[1] = {QkBitTerm_Z};
        uint32_t indices[1] = {i - 1};
        QkSparseTerm term_alpha = {coeff_half, 1, bit_terms, indices,
                                   num_qubits};
        qk_obs_add_term(terms, &term_alpha);
        // beta-spin: -0.5 * coeff * Z_{i + N}
        indices[0] = i - 1 + num_orbs;
        QkSparseTerm term_beta = {coeff_half, 1, bit_terms, indices,
                                  num_qubits};
        qk_obs_add_term(terms, &term_beta);

        obs = qk_obs_add(obs, terms);
      }
      else
      {
        QkSparseObservable *terms = qk_obs_zero(num_qubits);

        // we always add XX and YY with coeff_half
        double complex coeff_half = 0.5 * coeff;
        // alpha-spin
        uint32_t dist_ia = a - i + 1;
        uint32_t indices[dist_ia];
        QkBitTerm bit_terms[dist_ia];
        for (uint32_t k = 1; k < dist_ia - 1; k++)
        {
          bit_terms[k] = QkBitTerm_Z;
          indices[k] = i - 1 + k;
        }
        indices[0] = i - 1;
        indices[dist_ia - 1] = a - 1;
        bit_terms[0] = QkBitTerm_X;
        bit_terms[dist_ia - 1] = QkBitTerm_X;
        QkSparseTerm term_alpha_xx = {coeff_half, dist_ia, bit_terms, indices,
                                      num_qubits};
        qk_obs_add_term(terms, &term_alpha_xx);
        bit_terms[0] = QkBitTerm_Y;
        bit_terms[dist_ia - 1] = QkBitTerm_Y;
        QkSparseTerm term_alpha_yy = {coeff_half, dist_ia, bit_terms, indices,
                                      num_qubits};
        qk_obs_add_term(terms, &term_alpha_yy);
        // beta-spin
        for (uint32_t k = 1; k < dist_ia - 1; k++)
        {
          indices[k] = i - 1 + k + num_orbs;
        }
        indices[0] = i - 1 + num_orbs;
        indices[dist_ia - 1] = a - 1 + num_orbs;
        QkSparseTerm term_beta_yy = {coeff_half, dist_ia, bit_terms, indices,
                                     num_qubits};
        qk_obs_add_term(terms, &term_beta_yy);
        bit_terms[0] = QkBitTerm_X;
        bit_terms[dist_ia - 1] = QkBitTerm_X;
        QkSparseTerm term_beta_xx = {coeff_half, dist_ia, bit_terms, indices,
                                     num_qubits};
        qk_obs_add_term(terms, &term_beta_xx);

        obs = qk_obs_add(obs, terms);
      }
    }
    else
    {
      // FIXME: handle intermittent Z terms

      // 2-body term
      assert(i != 0);
      assert(a != 0);
      assert(j != 0);
      assert(b != 0);

      if (i == a && j == b)
      {
        if (i == j)
        {
          QkSparseObservable *terms = qk_obs_identity(num_qubits);

          double complex coeff_quart = 0.25 * coeff;
          double complex coeff_quart_minus = -0.25 * coeff;

          // 0.25 * coeff * Id
          terms = qk_obs_multiply(terms, &coeff_quart);

          // alpha-spin: -0.25 * coeff * Z_i
          QkBitTerm bit_terms[1] = {QkBitTerm_Z};
          uint32_t indices[1] = {i - 1};
          QkSparseTerm term_alpha = {coeff_quart_minus, 1, bit_terms, indices,
                                     num_qubits};
          qk_obs_add_term(terms, &term_alpha);

          // beta-spin: -0.25 * coeff * Z_{i + N}
          indices[0] = i - 1 + num_orbs;
          QkSparseTerm term_beta = {coeff_quart_minus, 1, bit_terms, indices,
                                    num_qubits};
          qk_obs_add_term(terms, &term_beta);

          // mixed-spin: 0.25 * coeff * Z_i @ Z_{i + N}
          QkBitTerm bit_terms_mixed[2] = {QkBitTerm_Z, QkBitTerm_Z};
          uint32_t indices_mixed[2] = {i - 1, i - 1 + num_orbs};
          QkSparseTerm term_mixed = {coeff_quart, 2, bit_terms_mixed,
                                     indices_mixed, num_qubits};
          qk_obs_add_term(terms, &term_mixed);

          obs = qk_obs_add(obs, terms);
        }
        else
        {
          QkSparseObservable *terms = qk_obs_identity(num_qubits);

          // coeff * Id
          terms = qk_obs_multiply(terms, &coeff);

          // single-Z terms: -0.5 * coeff * Z_k
          // for all k in {i, j, i + N, j + N}
          double complex coeff_single_z = -0.5 * coeff;
          QkBitTerm bit_terms_single_z[1] = {QkBitTerm_Z};
          uint32_t indices_single_z[1] = {i - 1};
          QkSparseTerm term_single_z = {coeff_single_z, 1, bit_terms_single_z,
                                        indices_single_z, num_qubits};
          qk_obs_add_term(terms, &term_single_z);
          indices_single_z[0] = j - 1;
          qk_obs_add_term(terms, &term_single_z);
          indices_single_z[0] = i - 1 + num_orbs;
          qk_obs_add_term(terms, &term_single_z);
          indices_single_z[0] = j - 1 + num_orbs;
          qk_obs_add_term(terms, &term_single_z);

          // double-Z terms: 0.25 * coeff * Z_k @ Z_l
          // for all (k, l) in {(i, j), (i + N, j), (i, j + N), (i + N, j + N)}
          double complex coeff_double_z = 0.25 * coeff;
          QkBitTerm bit_terms_double_z[2] = {QkBitTerm_Z, QkBitTerm_Z};
          uint32_t indices_double_z[2] = {i - 1, j - 1};
          QkSparseTerm term_double_z = {coeff_double_z, 2, bit_terms_double_z,
                                        indices_double_z, num_qubits};
          qk_obs_add_term(terms, &term_double_z);
          indices_double_z[0] = i - 1 + num_orbs;
          qk_obs_add_term(terms, &term_double_z);
          indices_double_z[1] = j - 1 + num_orbs;
          qk_obs_add_term(terms, &term_double_z);
          indices_double_z[0] = i - 1;
          qk_obs_add_term(terms, &term_double_z);

          obs = qk_obs_add(obs, terms);
        }
      }
      else if (i != a && j != b && i == j && a == b)
      {
        QkSparseObservable *terms = qk_obs_identity(num_qubits);

        // -0.5 * coeff * Id
        double complex coeff_half = -0.5 * coeff;
        terms = qk_obs_multiply(terms, &coeff_half);

        // single-Z terms: 0.25 * coeff * Z_k
        // for all k in {i, a, i + N, a + N}
        double complex coeff_single_z = 0.25 * coeff;
        QkBitTerm bit_terms_single_z[1] = {QkBitTerm_Z};
        uint32_t indices_single_z[1] = {i - 1};
        QkSparseTerm term_single_z = {coeff_single_z, 1, bit_terms_single_z,
                                      indices_single_z, num_qubits};
        qk_obs_add_term(terms, &term_single_z);
        indices_single_z[0] = a - 1;
        qk_obs_add_term(terms, &term_single_z);
        indices_single_z[0] = i - 1 + num_orbs;
        qk_obs_add_term(terms, &term_single_z);
        indices_single_z[0] = a - 1 + num_orbs;
        qk_obs_add_term(terms, &term_single_z);

        // double-Z terms: -0.25 * coeff * Z_k @ Z_l
        // for all (k, l) in {(i, a), (i + N, a + N)}
        double complex coeff_double_z = -0.25 * coeff;
        QkBitTerm bit_terms_double_z[2] = {QkBitTerm_Z, QkBitTerm_Z};
        uint32_t indices_double_z[2] = {i - 1, a - 1};
        QkSparseTerm term_double_z = {coeff_double_z, 2, bit_terms_double_z,
                                      indices_double_z, num_qubits};
        qk_obs_add_term(terms, &term_double_z);
        indices_double_z[0] = i - 1 + num_orbs;
        indices_double_z[1] = a - 1 + num_orbs;
        qk_obs_add_term(terms, &term_double_z);

        // mixed-spin terms: pairings of XX and YY on (i, a) etc.
        double complex coeff_mixed = 0.25 * coeff;
        QkBitTerm bit_terms_mixed[4] = {QkBitTerm_X, QkBitTerm_X, QkBitTerm_X,
                                        QkBitTerm_X};
        uint32_t indices_mixed[4] = {i - 1, a - 1, j - 1 + num_orbs,
                                     b - 1 + num_orbs};
        QkSparseTerm term_mixed = {coeff_mixed, 4, bit_terms_mixed,
                                   indices_mixed, num_qubits};
        qk_obs_add_term(terms, &term_mixed);
        bit_terms_mixed[2] = QkBitTerm_Y;
        bit_terms_mixed[3] = QkBitTerm_Y;
        qk_obs_add_term(terms, &term_mixed);
        bit_terms_mixed[0] = QkBitTerm_Y;
        bit_terms_mixed[1] = QkBitTerm_Y;
        qk_obs_add_term(terms, &term_mixed);
        bit_terms_mixed[2] = QkBitTerm_X;
        bit_terms_mixed[3] = QkBitTerm_X;
        qk_obs_add_term(terms, &term_mixed);

        obs = qk_obs_add(obs, terms);
      }

      // TODO: complete me!
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
