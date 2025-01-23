#include <assert.h>
#include <complex.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "qiskit.h"

int main() {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen("h2.fcidump", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    uintmax_t num_orbs = 2;
    uintmax_t num_qubits = 2 * num_orbs;

    SparseObservable* obs = qk_obs_zero(num_qubits);

    int line_co = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        line_co++;
        if (line_co < 5) {
            continue;
        }
        char* end = NULL;
        double c = strtod(line, &end);
        double complex coeff = c + 0.0 * I;
        uintmax_t i = strtoumax(end, &end, 10);
        uintmax_t a = strtoumax(end, &end, 10);
        uintmax_t j = strtoumax(end, &end, 10);
        uintmax_t b = strtoumax(end, &end, 10);

        if (i == 0) {
            // constant energy offset
            assert(a == 0);
            assert(j == 0);
            assert(b == 0);
            // coeff * Id
            BitTerm bit_terms[] = {};
            uint32_t indices[] = {};
            SparseTerm term = {&coeff, 0, bit_terms, indices, num_qubits};
            qk_obs_add_term(obs, &term);

        } else if (j == 0) {
            // 1-body term
            assert(b == 0);
            assert(i != 0);
            assert(a != 0);

            if (i == a) {
                SparseObservable* terms = qk_obs_identity(num_qubits);
                // coeff * Id
                terms = qk_obs_multiply(terms, &coeff);
                // alpha-spin: -0.5 * coeff * Z_i
                double complex coeff_half = -0.5 * coeff;
                BitTerm bit_terms[1] = {BitTerm_Z};
                uint32_t indices[1] = {i - 1};
                SparseTerm term_alpha = {&coeff_half, 1, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_alpha);
                // beta-spin: -0.5 * coeff * Z_{i + N}
                indices[0] = i - 1 + num_orbs;
                SparseTerm term_beta = {&coeff_half, 1, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_beta);

                obs = qk_obs_add(obs, terms);

            } else {
                SparseObservable* terms = qk_obs_zero(num_qubits);

                // we always add XX and YY with coeff_half
                double complex coeff_half = 0.5 * coeff;
                // alpha-spin
                uint32_t dist_ia = a - i + 1;
                uint32_t indices[dist_ia];
                BitTerm bit_terms[dist_ia];
                for (uint32_t k = 1; k < dist_ia - 1; k++) {
                    bit_terms[k] = BitTerm_Z;
                    indices[k] = i - 1 + k;
                }
                indices[0] = i - 1;
                indices[dist_ia - 1] = a - 1;
                bit_terms[0] = BitTerm_X;
                bit_terms[dist_ia - 1] = BitTerm_X;
                SparseTerm term_alpha_xx = {&coeff_half, dist_ia, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_alpha_xx);
                bit_terms[0] = BitTerm_Y;
                bit_terms[dist_ia - 1] = BitTerm_Y;
                SparseTerm term_alpha_yy = {&coeff_half, dist_ia, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_alpha_yy);
                // beta-spin
                for (uint32_t k = 1; k < dist_ia - 1; k++) {
                    indices[k] = i - 1 + k + num_orbs;
                }
                indices[0] = i - 1 + num_orbs;
                indices[dist_ia - 1] = a - 1 + num_orbs;
                SparseTerm term_beta_yy = {&coeff_half, dist_ia, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_beta_yy);
                bit_terms[0] = BitTerm_X;
                bit_terms[dist_ia - 1] = BitTerm_X;
                SparseTerm term_beta_xx = {&coeff_half, dist_ia, bit_terms, indices, num_qubits};
                qk_obs_add_term(terms, &term_beta_xx);

                obs = qk_obs_add(obs, terms);
            }

        } else {
            // FIXME: handle intermittent Z terms

            // 2-body term
            assert(i != 0);
            assert(a != 0);
            assert(j != 0);
            assert(b != 0);

            if (i == a && j == b) {
                if (i == j) {
                    SparseObservable* terms = qk_obs_identity(num_qubits);

                    double complex coeff_quart = 0.25 * coeff;
                    double complex coeff_quart_minus = -0.25 * coeff;

                    // 0.25 * coeff * Id
                    terms = qk_obs_multiply(terms, &coeff_quart);

                    // alpha-spin: -0.25 * coeff * Z_i
                    BitTerm bit_terms[1] = {BitTerm_Z};
                    uint32_t indices[1] = {i - 1};
                    SparseTerm term_alpha = {&coeff_quart_minus, 1, bit_terms, indices, num_qubits};
                    qk_obs_add_term(terms, &term_alpha);

                    // beta-spin: -0.25 * coeff * Z_{i + N}
                    indices[0] = i - 1 + num_orbs;
                    SparseTerm term_beta = {&coeff_quart_minus, 1, bit_terms, indices, num_qubits};
                    qk_obs_add_term(terms, &term_beta);

                    // mixed-spin: 0.25 * coeff * Z_i @ Z_{i + N}
                    BitTerm bit_terms_mixed[2] = {BitTerm_Z, BitTerm_Z};
                    uint32_t indices_mixed[2] = {i - 1, i - 1 + num_orbs};
                    SparseTerm term_mixed = {&coeff_quart, 2, bit_terms_mixed, indices_mixed, num_qubits};
                    qk_obs_add_term(terms, &term_mixed);

                    obs = qk_obs_add(obs, terms);

                } else {
                    SparseObservable* terms = qk_obs_identity(num_qubits);

                    // coeff * Id
                    terms = qk_obs_multiply(terms, &coeff);

                    // single-Z terms: -0.5 * coeff * Z_k
                    // for all k in {i, j, i + N, j + N}
                    double complex coeff_single_z = -0.5 * coeff;
                    BitTerm bit_terms_single_z[1] = {BitTerm_Z};
                    uint32_t indices_single_z[1] = {i - 1};
                    SparseTerm term_single_z = {&coeff_single_z, 1, bit_terms_single_z, indices_single_z, num_qubits};
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
                    BitTerm bit_terms_double_z[2] = {BitTerm_Z, BitTerm_Z};
                    uint32_t indices_double_z[2] = {i - 1, j - 1};
                    SparseTerm term_double_z = {&coeff_double_z, 2, bit_terms_double_z, indices_double_z, num_qubits};
                    qk_obs_add_term(terms, &term_double_z);
                    indices_double_z[0] = i - 1 + num_orbs;
                    qk_obs_add_term(terms, &term_double_z);
                    indices_double_z[1] = j - 1 + num_orbs;
                    qk_obs_add_term(terms, &term_double_z);
                    indices_double_z[0] = i - 1;
                    qk_obs_add_term(terms, &term_double_z);

                    obs = qk_obs_add(obs, terms);
                }

            } else if (i != a  && j != b && i == j && a == b) {
                    SparseObservable* terms = qk_obs_identity(num_qubits);

                    // -0.5 * coeff * Id
                    double complex coeff_half = -0.5 * coeff;
                    terms = qk_obs_multiply(terms, &coeff_half);

                    // single-Z terms: 0.25 * coeff * Z_k
                    // for all k in {i, a, i + N, a + N}
                    double complex coeff_single_z = 0.25 * coeff;
                    BitTerm bit_terms_single_z[1] = {BitTerm_Z};
                    uint32_t indices_single_z[1] = {i - 1};
                    SparseTerm term_single_z = {&coeff_single_z, 1, bit_terms_single_z, indices_single_z, num_qubits};
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
                    BitTerm bit_terms_double_z[2] = {BitTerm_Z, BitTerm_Z};
                    uint32_t indices_double_z[2] = {i - 1, a - 1};
                    SparseTerm term_double_z = {&coeff_double_z, 2, bit_terms_double_z, indices_double_z, num_qubits};
                    qk_obs_add_term(terms, &term_double_z);
                    indices_double_z[0] = i - 1 + num_orbs;
                    indices_double_z[1] = a - 1 + num_orbs;
                    qk_obs_add_term(terms, &term_double_z);

                    // mixed-spin terms: pairings of XX and YY on (i, a) etc.
                    double complex coeff_mixed = 0.25 * coeff;
                    BitTerm bit_terms_mixed[4] = {BitTerm_X, BitTerm_X, BitTerm_X, BitTerm_X};
                    uint32_t indices_mixed[4] = {i - 1, a - 1, j - 1 + num_orbs, b - 1 + num_orbs};
                    SparseTerm term_mixed = {&coeff_mixed, 4, bit_terms_mixed, indices_mixed, num_qubits};
                    qk_obs_add_term(terms, &term_mixed);
                    bit_terms_mixed[2] = BitTerm_Y;
                    bit_terms_mixed[3] = BitTerm_Y;
                    qk_obs_add_term(terms, &term_mixed);
                    bit_terms_mixed[0] = BitTerm_Y;
                    bit_terms_mixed[1] = BitTerm_Y;
                    qk_obs_add_term(terms, &term_mixed);
                    bit_terms_mixed[2] = BitTerm_X;
                    bit_terms_mixed[3] = BitTerm_X;
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
    qk_obs_print(obs);

    exit(EXIT_SUCCESS);
}
