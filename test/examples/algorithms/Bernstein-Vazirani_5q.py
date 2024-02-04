import pennylane as qml
import numpy as np

# dev = qml.device('default.qubit', wires=5)
dev = qml.device("qiskit.aer", wires=5)


@qml.qnode(dev)
def circuit():
    qml.RY(0.5 * np.pi, wires=0)
    qml.RY(0.5 * np.pi, wires=1)
    qml.RY(0.5 * np.pi, wires=2)
    qml.RY(0.5 * np.pi, wires=3)
    qml.RY(-0.5 * np.pi, wires=4)
    qml.CNOT(wires=[0, 4])
    qml.CNOT(wires=[1, 4])
    qml.CNOT(wires=[2, 4])
    qml.CNOT(wires=[3, 4])
    qml.RY(0.5 * np.pi, wires=0)
    qml.RY(0.5 * np.pi, wires=1)
    qml.RY(0.5 * np.pi, wires=2)
    qml.RY(0.5 * np.pi, wires=3)
    qml.RY(0.5 * np.pi, wires=4)
    # return qml.expval(qml.PauliZ(0)@qml.PauliZ(1)@qml.PauliZ(2)@qml.PauliZ(3))
    return qml.counts()


print(circuit())
print(dev._circuit.qasm(formatted=True))
