import os
import sys
from pprint import pprint
from datetime import datetime
import simplejson as json

# mapper directory
sys.path.append("../library")
import qcmapper

# 
import layout_generator
import util


if __name__ == "__main__":
    list_kisti_algorithms = [
                            "Bernstein-Vazirani_5q.qasm", 
                            # "Bernstein-Vazirani_5q_2.qasm",
                            # "CHSH1.qasm", "CHSH2.qasm", "CHSH3.qasm", 
                            # "CHSH4.qasm",
                            "ae_indep_qiskit_5.qasm",
                            # "dj_indep_qiskit_5.qasm",
                            # "ghz_indep_qiskit_5.qasm",
                            # "graphstate_indep_qiskit_5.qasm",
                            # "grover-noancilla_indep_qiskit_5.qasm",
                            # "grover-v-chain_indep_qiskit_5.qasm",
                            # "portfolioqaoa_indep_qiskit_5.qasm",
                            # "portfoliovqe_indep_qiskit_5.qasm",
                            # "qaoa_indep_qiskit_5.qasm",
                            "qft_indep_qiskit_5.qasm",
                            # "qftentangled_indep_qiskit_5.qasm",
                            # "qgan_indep_qiskit_5.qasm",
                            # "qpeexact_indep_qiskit_5.qasm",
                            # "qpeinexact_indep_qiskit_5.qasm",
                            # "qwalk-noancilla_indep_qiskit_5.qasm",
                            # "qwalk-v-chain_indep_qiskit_5.qasm",
                            # "realamprandom_indep_qiskit_5.qasm",
                            # "su2random_indep_qiskit_5.qasm",
                            # "twolocalrandom_indep_qiskit_5.qasm",
                            # "vqe_indep_qiskit_5.qasm",
                            # "wstate_indep_qiskit_5.qasm",
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
                      "calibration" : True,
                      "iteration": 10, 
                      "optimal_criterion": "fidelity",
                      "initial_mapping_option": "random", 
                      "mapper": "SABRE",
                      "cost": "lap",
                     }

    ##########################
    # instruction format options
    # 
    #   1. gate_angle : True -> ['U(theta, phi, lambda)', 'qubit']
    #                   False -> ['U', 'theta', 'phi', 'lambda', 'qubit']
    
    instruction_format = {"gate_angle": True}

    ##########################
    # quantum chip information

    # 1. user quantum chip (provided by the file path)
    path_qchip = os.path.join("examples/quantum_chips", "ibmq_16_melbourne.json")
    user_chip = True
    
    # 2. artificially generated quantum chip by using a delivered package (layout_generator)
    # architecture : {0: all-to-all, 
    #                 2: 2D rectangular, 
    #                 23: 2D rectangular but having triangle face
    #                 3: 3D rectangular }
    
    # import layout_generator
    # chip_dimension = {"height": 3, "width": 3, "length": 1}
    # qchip = layout_generator.generate_regular_qchip_architecture("examples/quantum_chips", chip_dimension, 
    #         architecture=2)
    # user_chip = False
    # path_qchip = qchip.get("result_file")

    dir_jobs = os.path.join("examples", "jobs")
    if not os.path.exists(dir_jobs):
        os.makedirs(dir_jobs)

    now = datetime.now()

    for algorithm in list_kisti_algorithms:
        print(algorithm)
        path_qasm = os.path.join("examples/algorithms", algorithm)
        ret = qcmapper.map_circuit(path_qasm, path_qchip, option=synthesis_option,
                                                          format=instruction_format)
        pprint(ret)

        # writing a circuit into a file
        algorithm_name = os.path.splitext(algorithm)[0]
        qchip_name = os.path.splitext(os.path.basename(path_qchip))[0]
        current_time = "{}-{}-{}".format(now.date(),now.hour,now.minute)
        task_id = "{}-{}-{}".format(algorithm_name, qchip_name, current_time)
        
        dir_task = os.path.join(dir_jobs, task_id)
        os.makedirs(dir_task)
        file_circuit = os.path.join(dir_task, "circuit.json")
        
        with open(file_circuit, "w") as outfile:
            json.dump(ret, outfile, sort_keys=True, indent=4, separators=(',', ':'))

        # functions to display the qubit movings
        # to see the position of qubit during the circuit
        util.display_qubit_movements(ret.get("system_code"), qchip=ret.get("qchip"))
        