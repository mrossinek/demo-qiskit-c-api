import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import excitation_preserving
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparseObservable, SparsePauliOp

import cextension

if __name__ == "__main__":
    obs = cextension.qubit_observable()
    print("-- SparseObservable")
    print(obs)

    circuit = QuantumCircuit(obs.num_qubits)
    circuit.x(0)
    circuit.compose(excitation_preserving(obs.num_qubits), inplace=True)
    values = 2 * np.pi * np.random.random((10, circuit.num_parameters))

    # primitives V2 don't support SparseObservable yet
    spo = SparsePauliOp.from_sparse_observable(obs)
    pub = (circuit, spo, values)
    result = StatevectorEstimator().run([pub]).result()[0]

    print("-- Estimator results")
    print(result.data.evs)
