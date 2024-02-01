# How to use the sample code *qcmapper\_test.py*

## How to execute
```
python qcmapper_test.py
```
```
import os
import sys

# package import
import library.qcmapper as qcmapper
from pprint import pprint

if __name__ == "__main__":
    list_kisti_algorithms = [
                            "Bernstein-Vazirani_5q.qasm", 
                            "Bernstein-Vazirani_5q_2.qasm",
                            "CHSH1.qasm",
                            "CHSH2.qasm",
                            "CHSH3.qasm", 
                            "CHSH4.qasm"
                            ]
    # mapping options
    synthesis_option={"allow_swap": True, "random_seed": 0,
                      "calibration" : True,
                      "iteration": 10, "cost": "lap", "optimal_criterion": "number_gates",
                      "initial_mapping_option": "random", "mapper": "SABRE"}

    # example quantum chip
    path_qchip = os.path.join("examples/quantum_chips", "ibmq_16_melbourne.json")

    for algorithm in list_kisti_algorithms:
    	path_qasm = os.path.join("examples/algorithms/KISTI_sample", algorithm)

      # call the mapper 
    	ret = qcmapper.map_circuit(path_qasm, path_qchip, option=synthesis_option)

    	pprint(ret)
```
