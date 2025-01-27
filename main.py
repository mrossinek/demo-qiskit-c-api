import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import excitation_preserving
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparseObservable, SparsePauliOp

import cmod

def to_sparse_pauli_op(obs: SparseObservable) -> SparsePauliOp:
    bitlabels = {1: "Z", 2: "X", 3: "Y"}
    to_bitstring = lambda bit_terms: "".join(bitlabels[b] for b in bit_terms)

    sparse_labels = [(to_bitstring(term.bit_terms), term.indices.astype(int).tolist(), term.coeff) for term in obs]
    print(sparse_labels)
    print(obs.num_qubits)
    return SparsePauliOp.from_sparse_list(sparse_labels, obs.num_qubits)

if __name__ == "__main__":
    capsule = cmod.get_qubit_observable()
    obs = SparseObservable.from_pycapsule(capsule)
    print("-- SparseObservable")
    print(obs)

    spo = to_sparse_pauli_op(obs)
    print("-- SparsePauliOp")
    print(spo)

    circuit = QuantumCircuit(obs.num_qubits)
    circuit.x(0)
    circuit.compose(excitation_preserving(spo.num_qubits), inplace=True)
    values = 2 * np.pi * np.random.random((10, circuit.num_parameters))

    pub = (circuit, spo, values)
    result = StatevectorEstimator().run([pub]).result()[0]

    print("-- Estimator results")
    print(result.data.evs)
