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
            // 2-body term
            assert(i != 0);
            assert(a != 0);
            assert(j != 0);
            assert(b != 0);

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
