// for Python interface
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// main C script
#include <assert.h>
#include <complex.h>
#include <inttypes.h>
#include <qiskit.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

QkSparseObservable *jw_term(uintmax_t num_qubits, uintmax_t site, bool creation) {
    QkSparseObservable *op = qk_obs_zero(num_qubits);

    uint32_t indices[site];
    QkBitTerm x_bits[site];
    QkBitTerm y_bits[site]; // could just allocate 1 of these

    for (uintmax_t i = 0; i < site - 1; i++) {
        indices[i] = i;
        x_bits[i] = QkBitTerm_Z;
        y_bits[i] = QkBitTerm_Z;
    }
    indices[site - 1] = site - 1;
    x_bits[site - 1] = QkBitTerm_X;
    y_bits[site - 1] = QkBitTerm_Y;

    complex double x_coeff = 0.5;
    complex double y_coeff = creation ? -0.5 * I : 0.5 * I;

    QkSparseTerm x = {x_coeff, site, x_bits, indices, num_qubits};
    QkSparseTerm y = {y_coeff, site, y_bits, indices, num_qubits};

    qk_obs_add_term(op, &x);
    qk_obs_add_term(op, &y);

    return op;
}

void add_one_body(QkSparseObservable *obs, uintmax_t num_qubits, uintmax_t num_orbs,
                  complex double coeff, uintmax_t i, uintmax_t a) {
    QkSparseObservable *a_i = jw_term(num_qubits, i, true);
    QkSparseObservable *a_a = jw_term(num_qubits, a, false);
    QkSparseObservable *out = qk_obs_compose(a_a, a_i);
    qk_obs_multiply_inplace(out, &coeff);

    qk_obs_append(obs, out);

    qk_obs_free(a_i);
    qk_obs_free(a_a);
    qk_obs_free(out);
}

void add_one_body_spin(QkSparseObservable *obs, uintmax_t num_qubits, uintmax_t num_orbs,
                       complex double coeff_a, complex double coeff_b, uintmax_t i, uintmax_t a) {
    add_one_body(obs, num_qubits, num_orbs, coeff_a, i, a);
    add_one_body(obs, num_qubits, num_orbs, coeff_b, i + num_orbs, a + num_orbs);
    if (a != i) {
        add_one_body(obs, num_qubits, num_orbs, coeff_a, a, i);
        add_one_body(obs, num_qubits, num_orbs, coeff_b, a + num_orbs, i + num_orbs);
    }
}

void add_two_body(QkSparseObservable *obs, uintmax_t num_qubits, uintmax_t num_orbs,
                  complex double coeff, uintmax_t i, uintmax_t a, uintmax_t j, uintmax_t b) {
    QkSparseObservable *a_i = jw_term(num_qubits, i, true);
    QkSparseObservable *a_a = jw_term(num_qubits, a, false);

    QkSparseObservable *a_j = jw_term(num_qubits, j, true);
    QkSparseObservable *a_b = jw_term(num_qubits, b, false);

    // math world: create(i) create(j) ann(b) ann(a)
    QkSparseObservable *right = qk_obs_compose(a_a, a_b);
    QkSparseObservable *left = qk_obs_compose(a_j, a_i);
    QkSparseObservable *out = qk_obs_compose(right, left);

    qk_obs_multiply_inplace(out, &coeff);

    qk_obs_append(obs, out);

    qk_obs_free(a_i);
    qk_obs_free(a_j);
    qk_obs_free(a_a);
    qk_obs_free(a_b);
    qk_obs_free(out);
    qk_obs_free(right);
    qk_obs_free(left);
}

void add_two_body_permute(QkSparseObservable *obs, uintmax_t num_qubits, uintmax_t num_orbs,
                          complex double coeff, bool swap_ia_jb, uintmax_t i, uintmax_t a,
                          uintmax_t j, uintmax_t b) {
    add_two_body(obs, num_qubits, num_orbs, coeff, i, a, j, b);
    if (b > j)
        add_two_body(obs, num_qubits, num_orbs, coeff, i, a, b, j);
    if (a > i) {
        add_two_body(obs, num_qubits, num_orbs, coeff, a, i, j, b);
        if (b > j)
            add_two_body(obs, num_qubits, num_orbs, coeff, a, i, b, j);
    }

    if (swap_ia_jb) {
        // swap i with j and a with b
        add_two_body(obs, num_qubits, num_orbs, coeff, j, b, i, a);
        if (a > i)
            add_two_body(obs, num_qubits, num_orbs, coeff, j, b, a, i);
        if (b > j) {
            add_two_body(obs, num_qubits, num_orbs, coeff, b, j, i, a);
            if (a > i)
                add_two_body(obs, num_qubits, num_orbs, coeff, b, j, a, i);
        }
    }
}

void add_two_body_spin_pure(QkSparseObservable *obs, uintmax_t num_qubits, uintmax_t num_orbs,
                            complex double coeff_a, complex double coeff_b, bool swap_ia_jb,
                            uintmax_t i, uintmax_t a, uintmax_t j, uintmax_t b) {
    add_two_body_permute(obs, num_qubits, num_orbs, coeff_a, swap_ia_jb, i, a, j, b);
    add_two_body_permute(obs, num_qubits, num_orbs, coeff_b, swap_ia_jb, i + num_orbs, a + num_orbs,
                         j + num_orbs, b + num_orbs);
}

void _inflate_index(uintmax_t ia, uintmax_t size, uintmax_t *i, uintmax_t *a) {
    // reverse an index of a flattened upper-triangular array (ia) to its original index pair (i, a)
    *i = 0;
    *a = 0;
    for (uintmax_t n = 0; n < size; n++) {
        uintmax_t pq = (n + 1) * (n + 2) / 2;
        if (pq > ia)
            break;
        *i += 1;
    }
    *a = ia - *i * (*i + 1) / 2;
}

uintmax_t _deflate_index(uintmax_t i, uintmax_t a) {
    // computes the flattened upper-triangular array index (ia) from an index pair (i, a)
    if (a < i) {
        uintmax_t tmp = a;
        a = i;
        i = tmp;
    }
    return a * (a - 1) / 2 + i - 1;
}

QkSparseObservable *get_molecular_hamiltonian(char *filename) {
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(filename, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    uintmax_t num_orbs = 0;
    uintmax_t num_qubits = 0;
    uintmax_t num_pairs = 0;

    double *one_body;
    double *two_body;
    double *one_body_b;
    double *two_body_bb;
    double *two_body_ab;
    bool unrestricted_spin = false;

    QkSparseObservable *obs;

    bool finished_header = false;
    while ((read = getline(&line, &len, fp)) != -1) {
        if (!finished_header) {
            if (num_orbs == 0) {
                // we have not detected the number of orbitals, yet
                // FIXME: this assumes that NORB is the first parameter in the
                // namelist, but this is not guaranteed to be the case.
                sscanf(line, "&FCI NORB = %ju,", &num_orbs);
                num_qubits = 2 * num_orbs;
                obs = qk_obs_zero(num_qubits);
                // we know that the header is not finished yet, so we can skip
                // the next check
                continue;
            }
            double c;
            unsigned int num_chars;
            // once the line starts with a float, the header is finished
            sscanf(line, "%le%n", &c, &num_chars);
            if (num_chars < 100) {
                finished_header = true;
                // allocate memory for one- and two-body coefficients arrays
                num_pairs = num_orbs * (num_orbs + 1) / 2;
                // NOTE: we use calloc to ensure that our arrays are initialized with zeros
                one_body = calloc(num_pairs, sizeof(double));
                // The 2-body terms are 8-fold symmetric, meaning we only need to store the upper-
                // triangular of the upper-triangular of the 4D matrix (hence this array size).
                two_body = calloc(num_pairs * (num_pairs + 1) / 2, sizeof(double));
            } else {
                continue;
            }
        }

        char *end = NULL;
        double c = strtod(line, &end);
        double complex coeff = c + 0.0 * I;
        uintmax_t i = strtoumax(end, &end, 10);
        uintmax_t a = strtoumax(end, &end, 10);
        uintmax_t j = strtoumax(end, &end, 10);
        uintmax_t b = strtoumax(end, &end, 10);
        uintmax_t ia;
        uintmax_t jb;
        uintmax_t iajb;

        bool ia_beta = false;
        bool jb_beta = false;
        double *hij = one_body;
        double *hijkl = two_body;

        if (i == 0 && j == 0 && a == 0 && b == 0) {
            uint32_t inds[] = {};
            QkBitTerm bits[] = {};
            QkSparseTerm energy = {coeff, 0, bits, inds, num_qubits};
            qk_obs_add_term(obs, &energy);

        } else if (j == 0 && b == 0) {
            ia_beta = i > num_orbs || a > num_orbs;
            if (ia_beta) {
                if (!unrestricted_spin) {
                    unrestricted_spin = true;
                    // NOTE: we use calloc to ensure that our arrays are initialized with zeros
                    one_body_b = calloc(num_pairs, sizeof(double));
                    two_body_bb = calloc(num_pairs * (num_pairs + 1) / 2, sizeof(double));
                    // NOTE: the mixed-spin integrals are only 4-fold symmetric, hence the different
                    // size of the array
                    two_body_ab = calloc(num_pairs * num_pairs, sizeof(double));
                }

                i -= num_orbs;
                a -= num_orbs;
                hij = one_body_b;
            }

            ia = _deflate_index(i, a);
            hij[ia] = c;

        } else {
            ia_beta = i > num_orbs || a > num_orbs;
            jb_beta = j > num_orbs || b > num_orbs;
            bool mixed_spin = false;
            if (ia_beta || jb_beta) {
                if (!unrestricted_spin) {
                    unrestricted_spin = true;
                    // NOTE: we use calloc to ensure that our arrays are initialized with zeros
                    one_body_b = calloc(num_pairs, sizeof(double));
                    two_body_bb = calloc(num_pairs * (num_pairs + 1) / 2, sizeof(double));
                    // NOTE: the mixed-spin integrals are only 4-fold symmetric, hence the different
                    // size of the array
                    two_body_ab = calloc(num_pairs * num_pairs, sizeof(double));
                }

                j -= num_orbs;
                b -= num_orbs;

                if (ia_beta && jb_beta) {
                    hijkl = two_body_bb;
                    i -= num_orbs;
                    a -= num_orbs;
                } else {
                    mixed_spin = true;
                    hijkl = two_body_ab;
                }
            }

            ia = _deflate_index(i, a);
            jb = _deflate_index(j, b);
            // NOTE: we must account for the 1-based indexing in _deflate_index
            iajb = mixed_spin ? ia * num_pairs + jb : _deflate_index(jb + 1, ia + 1);
            hijkl[iajb] = 0.5 * c;
        }
    }

    uintmax_t i[1];
    uintmax_t a[1];
    uintmax_t j[1];
    uintmax_t b[1];
    double complex coeff_a;
    double complex coeff_b;
    for (uintmax_t ia = 0; ia < num_pairs; ia++) {
        _inflate_index(ia, num_orbs, a, i);
        // we need to account for the 1-based indices
        *i += 1;
        *a += 1;

        coeff_a = one_body[ia] + 0.0 * I;
        coeff_b = unrestricted_spin ? one_body_b[ia] : coeff_a + 0.0 * I;
        add_one_body_spin(obs, num_qubits, num_orbs, coeff_a, coeff_b, *i, *a);

        for (uintmax_t jb = unrestricted_spin ? 0 : ia; jb < num_pairs; jb++) {
            _inflate_index(jb, num_orbs, b, j);
            // we need to account for the 1-based indices
            *j += 1;
            *b += 1;

            // NOTE: we only add the pure-spin 2-body terms when jb>=ia because only then does
            // the additonal symmetry hold
            if (jb >= ia) {
                // NOTE: we must account for the 1-based indexing in _deflate_index
                uintmax_t iajb = _deflate_index(jb + 1, ia + 1);
                coeff_a = two_body[iajb] + 0.0 * I;
                coeff_b = unrestricted_spin ? two_body_bb[iajb] : coeff_a + 0.0 * I;
                add_two_body_spin_pure(obs, num_qubits, num_orbs, coeff_a, coeff_b, (ia != jb), *i,
                                       *a, *j, *b);
            }

            if (unrestricted_spin) {
                uintmax_t iajb = ia * num_pairs + jb;
                double complex coeff = two_body_ab[iajb] + 0.0 * I;
                add_two_body_permute(obs, num_qubits, num_orbs, coeff, false, *i, *a, *j + num_orbs,
                                     *b + num_orbs);
                add_two_body_permute(obs, num_qubits, num_orbs, coeff, false, *j + num_orbs,
                                     *b + num_orbs, *i, *a);
            } else {
                // NOTE: coeff_a here is always well defined because this clause only gets executed
                // for jb>=ia due to the starting point's dependene on `unrestricted_spin`
                add_two_body_permute(obs, num_qubits, num_orbs, coeff_a, (ia != jb), *i + num_orbs,
                                     *a + num_orbs, *j, *b);
                add_two_body_permute(obs, num_qubits, num_orbs, coeff_a, (ia != jb), *i, *a,
                                     *j + num_orbs, *b + num_orbs);
            }
        }
    }

    fclose(fp);
    if (line)
        free(line);

    QkSparseObservable *canonicalized = qk_obs_canonicalize(obs, 1e-10);
    qk_obs_free(obs);

    return canonicalized;
}

void add_nn_interaction(QkSparseObservable *obs, double J, uint32_t i, uint32_t j,
                        uint32_t num_qubits) {

    QkBitTerm interaction[2] = {QkBitTerm_X, QkBitTerm_X};
    uint32_t indices[2];
    indices[0] = i * num_qubits + j;

    complex double coeff = J + 0.0 * I;
    QkSparseTerm term = {coeff, 2, interaction, indices, num_qubits * num_qubits};

    if (i > 0) {
        indices[1] = (i - 1) * num_qubits + j;
        qk_obs_add_term(obs, &term);
    }
    if (i < num_qubits - 1) {
        indices[1] = (i + 1) * num_qubits + j;
        qk_obs_add_term(obs, &term);
    }
    if (j > 0) {
        indices[1] = i * num_qubits + j - 1;
        qk_obs_add_term(obs, &term);
    }
    if (j < num_qubits - 1) {
        indices[1] = i * num_qubits + j + 1;
        qk_obs_add_term(obs, &term);
    }
}

QkSparseObservable *get_ising_lattice(uint32_t size) {
    uint32_t num_qubits = size * size;
    QkSparseObservable *obs = qk_obs_zero(num_qubits);
    double J = 0.5;
    double h = -1;

    QkBitTerm x[1] = {QkBitTerm_X};
    QkBitTerm zz[2] = {QkBitTerm_Z, QkBitTerm_Z};

    uint32_t field_idx[1];
    uint32_t inter_idx[2];
    QkSparseTerm field_term = {h, 1, x, field_idx, num_qubits};
    QkSparseTerm inter_term = {J, 2, zz, inter_idx, num_qubits};

    for (uint32_t i = 0; i < size; i++) {
        for (uint32_t j = 0; j < size; j++) {
            // set X term
            field_idx[0] = i * size + j;
            qk_obs_add_term(obs, &field_term);

            // set interaction terms
            inter_idx[0] = i * size + j;
            if (i > 0) {
                inter_idx[1] = (i - 1) * size + j;
                qk_obs_add_term(obs, &inter_term);
            }
            if (i < size - 1) {
                inter_idx[1] = (i + 1) * size + j;
                qk_obs_add_term(obs, &inter_term);
            }
            if (j > 0) {
                inter_idx[1] = i * size + j - 1;
                qk_obs_add_term(obs, &inter_term);
            }
            if (j < size - 1) {
                inter_idx[1] = i * size + j + 1;
                qk_obs_add_term(obs, &inter_term);
            }
        }
    }

    return obs;
}

static PyObject *cextension_molecular_hamiltonian(PyObject *self, PyObject *args) {
    char *filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    QkSparseObservable *obs = get_molecular_hamiltonian(filename);
    PyObject *py_obs = qk_obs_to_python(obs);
    return py_obs;
}

static PyObject *cextension_ising_observable(PyObject *self, PyObject *args) {
    unsigned int num_qubits;
    if (!PyArg_ParseTuple(args, "I", &num_qubits)) {
        return NULL;
    }

    QkSparseObservable *obs = get_ising_lattice((uint32_t)num_qubits);
    PyObject *py_obs = qk_obs_to_python(obs);
    return py_obs;
}

static PyMethodDef CExtMethods[] = {
    {"molecular_hamiltonian", cextension_molecular_hamiltonian, METH_VARARGS,
     "Get the qubit observable"},
    {"ising_observable", cextension_ising_observable, METH_VARARGS,
     "Get a Ising Hamiltonian on a lattice"},
    {NULL, NULL, 0, NULL}, // sentinel
};

static struct PyModuleDef cextension = {
    PyModuleDef_HEAD_INIT,
    "cextension", // module name
    NULL,         // docs
    -1,           // keep the module state in global variables
    CExtMethods,
};

PyMODINIT_FUNC PyInit_cextension(void) { return PyModule_Create(&cextension); }

int main(int argc, char *argv[]) {
    PyStatus status;
    PyConfig config;
    PyConfig_InitPythonConfig(&config);

    /* Add a built-in module, before Py_Initialize */
    if (PyImport_AppendInittab("cextension", PyInit_cextension) == -1) {
        fprintf(stderr, "Error: could not extend in-built modules table\n");
        exit(1);
    }

    /* Pass argv[0] to the Python interpreter */
    status = PyConfig_SetBytesString(&config, &config.program_name, argv[0]);
    if (PyStatus_Exception(status)) {
        goto exception;
    }

    /* Initialize the Python interpreter.  Required.
       If this step fails, it will be a fatal error. */
    status = Py_InitializeFromConfig(&config);
    if (PyStatus_Exception(status)) {
        goto exception;
    }
    PyConfig_Clear(&config);

    /* Optionally import the module; alternatively,
       import can be deferred until the embedded script
       imports it. */
    PyObject *pmodule = PyImport_ImportModule("cextension");
    if (!pmodule) {
        PyErr_Print();
        fprintf(stderr, "Error: could not import module 'cextension'\n");
    }

    return 0;

exception:
    PyConfig_Clear(&config);
    Py_ExitStatusException(status);
}
