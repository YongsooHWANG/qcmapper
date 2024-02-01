import pennylane as qml

dev = qml.device("qiskit.aer", wires=5)


@qml.qnode(dev)
def circuit1(x):
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    qml.RY(x, wires=0)
    return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))


@qml.qnode(dev)
def circuit2(x):
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    qml.RY(x, wires=0)
    qml.Hadamard(wires=0)
    return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))


@qml.qnode(dev)
def circuit3(x):
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    qml.RY(x, wires=0)
    qml.Hadamard(wires=1)
    return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))


@qml.qnode(dev)
def circuit4(x):
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    qml.RY(x, wires=0)
    qml.Hadamard(wires=0)
    qml.Hadamard(wires=1)
    return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))


def CHSH(x):
    a = circuit1(x) - circuit2(x) + circuit3(x) + circuit4(x)
    return a


import numpy as np

theta = [0.1 * i * np.pi for i in range(21)]

for i in theta:
    c = CHSH(i)
    print(c)


print(dev._circuit.qasm(formatted=True))
