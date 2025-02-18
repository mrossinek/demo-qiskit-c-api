from qiskit import transpile
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparseObservable
from qiskit.transpiler.passes import HLSConfig

import cextension

if __name__ == "__main__":
    obs = cextension.ising_observable(10)

    evo = PauliEvolutionGate(obs)
    circuit = QuantumCircuit(obs.num_qubits)
    circuit.append(evo, circuit.qubits)

    basis_gates = ["sx", "x", "rz", "cx"]
    isa = transpile(circuit, basis_gates=basis_gates)
    print("\nOptimization level 2")
    print("ops", isa.count_ops())
    print("2q depth", isa.depth(lambda inst: len(inst.qubits) > 1))

    config = HLSConfig(PauliEvolution=[("rustiq", {"preserve_order": False})])
    isa_rustiq = transpile(circuit, basis_gates=basis_gates, hls_config=config)
    print("\nwith Rustiq plugin")
    print("ops", isa_rustiq.count_ops())
    print("2q depth", isa_rustiq.depth(lambda inst: len(inst.qubits) > 1))
