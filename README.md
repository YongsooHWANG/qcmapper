# Quantum Circuit Mapper
- This project develops a quantum circuit mapping algorithm for computational algorithms to execute them on qubit-connectivity-constrained quantum computing HW.
- This project provides two functions: circuit mappings based only on qubit connectivity and the calibration data such as gate fidelity and time.
- For the connectivity resolving circuit mapping, we first implemented the [SABRE algorithm](https://dl.acm.org/doi/10.1145/3297858.3304023) and, second, developed a fast circuit mapping based on dijkstra's shortest path algorithm.
- For the calibration-aware mapping, we have redesigned and implemented the SABRE and the dijkstra based mapping.

## Environment
- Language :  Python3
- OS:  Ubuntu 20.04 

## Prerequisites
To run the project successfully, you need to install the following packages included in "requirements.txt" after installation.
- *simplejson*, *pandas*, *networkx*, *progress*, *icecream*, *qubitmapping*

```
pip install -r requirements.txt
```
Note that the packages *qubitmapping* is developed by Y.Hwang for this project.

## Installation
We encourage installing this project by cloning the source code from GitHub server.
But, we are working now that this project can be installed via pip. 
```
bash
git clone https://github.com/YongsooHWANG/qcmapper.git
```

## Usage

- The detailed usage and the options for the execution will be provided soon.
- For the sample demonstration, please see [Demo.md](docs/Demo.md).
- The test code and examples of algorithms and quantum chips are included in the directory **test**.
- This circuit mapper works well with the openqasm formatted quantum algorithm. Please see an example code [Bernstein-Vazirani_5q.qasm](test/examples/algorithms/Bernstein-Vazirani_5q.qasm).
- The quantum chip information file should be provided as a json format. Please see an example file [ibmq_16_melbourne.json](test/examples/quantum_chips/ibmq_16_melbourne.json). The template will be provided soon with a detailed explanation.

## Authors
- Yongsoo Hwang (ETRI, Quantum Technology Research Department), yhwang@etri.re.kr

## Reference
-

## License
This project is licensed under the [BSD-3-Clause](/docs/LICENSE.md)


