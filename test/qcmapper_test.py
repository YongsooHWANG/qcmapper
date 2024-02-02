import os
import sys
from icecream import ic
sys.path.append("../library/")

import qcmapper

from pprint import pprint


import layout_generator
import util
# from icecream import ic


if __name__ == "__main__":
    list_kisti_algorithms = [
                            "Bernstein-Vazirani_5q.qasm", 
                            # "Bernstein-Vazirani_5q_2.qasm",
                            # "CHSH1.qasm", "CHSH2.qasm", "CHSH3.qasm", 
                            # "CHSH4.qasm"
                            ]
    
    ##########################
    # circuit mapping options
    # 
    #   1. allow_swap : True -> inserting "SWAP" gate
    #                   False -> inserting 3 "CNOT" gates rather than 1 "SWAP" gate
    #   2. calibration : True -> circuit mapping based on calibration data (cnot gate time, fidelity,..)
    #                    False -> circuit mapping based on qubit connectivity only
    #   3. iteration : the rounds of circuit mappings (the more, the compact sized circuit)
    #   4. optimal_criterion : {"number_gates", "circuit_depth"} for non-calibration circuit mapping
    #                          {"fidelity", "time"} for calibration circuit mapping
    #   5. initial_mapping_option : "random" for nisq computational algorithms
    #   6. mapper : {"SABRE", "dijkstra"}
    #               SABRE : for more compact-sized circuit
    #               dijkstra : for faster circuit mapping
    #   7. cost : {"lap", "nnc"} for only SABRE mapper
    
    synthesis_option={"allow_swap": True, 
                      "calibration" : False,
                      "iteration": 2, 
                      "optimal_criterion": "fidelity",
                      "initial_mapping_option": "random", 
                      "mapper": "SABRE",
                      "cost": "nnc"
                     }

    ##########################
    # quantum chip information

    # 1. user quantum chip (provided by the file path)
    # path_qchip = os.path.join("examples/quantum_chips", "ibmq_16_melbourne.json")
    
    # 2. artificially generated quantum chip by using a delivered package (layout_generator)
    # architecture : {0: all-to-all, 
    #                 2: 2D rectangular, 
    #                 23: 2D rectangular but having triangle face
    #                 3: 3D rectangular }
    import layout_generator
    qchip = layout_generator.generate_regular_qchip_architecture("examples/quantum_chips", {"height": 3, "width": 3, "length": 2}, 
            architecture=2)
    path_qchip = qchip.get("result_file")

    for algorithm in list_kisti_algorithms:
        path_qasm = os.path.join("examples/algorithms", algorithm)
        ret = qcmapper.map_circuit(path_qasm, path_qchip, option=synthesis_option)

        pprint(ret)

        # function to display qubit movements on the chip
        chip_size = len(ret.get("qchip").get("qubit_connectivity"))
        util.display_qubit_movements(ret.get("system_code"), {"width": 1, "height": chip_size, "lenght": 1},)
        
