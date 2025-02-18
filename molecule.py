import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import excitation_preserving
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparseObservable, SparsePauliOp
from qiskit_nature.second_q.formats import fcidump_to_problem
from qiskit_nature.second_q.formats.fcidump import FCIDump
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

import cextension

if __name__ == "__main__":
    filename = "h2.fcidump"
    obs = cextension.molecular_hamiltonian(filename)
    print("-- SparseObservable")
    print(obs)

    # primitives V2 don't support SparseObservable yet...
    # ... and we want to compare against Qiskit Nature's output
    spo = SparsePauliOp.from_sparse_observable(obs)

    fcidump = FCIDump.from_file(filename)
    problem = fcidump_to_problem(fcidump)
    mapper = JordanWignerMapper()
    hamil = problem.hamiltonian.second_q_op()
    hamil += sum(problem.hamiltonian.constants.values()) * FermionicOp.one()
    op = mapper.map(hamil)

    assert spo.equiv(op)

    circuit = QuantumCircuit(obs.num_qubits)
    circuit.x(0)
    circuit.compose(excitation_preserving(obs.num_qubits), inplace=True)
    values = 2 * np.pi * np.random.random((10, circuit.num_parameters))

    pub = (circuit, spo, values)
    result = StatevectorEstimator().run([pub]).result()[0]

    print("-- Estimator results")
    print(result.data.evs)
