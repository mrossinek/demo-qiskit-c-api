import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import excitation_preserving
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparseObservable, SparsePauliOp

import cmod

if __name__ == "__main__":
    capsule = cmod.get_qubit_observable()
    obs = SparseObservable.from_pycapsule(capsule)
    print("-- SparseObservable")
    print(obs)

    spo = SparsePauliOp.from_sparse_observable(obs)
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
